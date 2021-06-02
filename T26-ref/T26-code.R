library(DescTools)

# Sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd() # is used in code - need to Source script, not Run
source("../utils.R")

# MAIN SCRIPT
#05/29/2021 T26 - Reference DUP240s
set.seed(54101)
vec <- sample(paste_expand(c("classII", "CLI1", "CLI2", "CLI3", "CLI4", "CLI5", "IFT1", "IFT2", "YJM1", "YJM2", "YJM3",
                               "YJM4", "YJM5", "YAR029W", "YPS1", "PWF1", "YAR023C", "YAR027W", "YAR028W",
                               "PRM9", "MST28", "YCR007C", "YHL044W", "YJM2+4"), "-", 1:4)
              , 96, replace=FALSE)
t(matrix(vec, nrow=8, dimnames = list(c(paste0("Column ", 3:10)), c(paste_expand("Plate ", 1:2, " ", c("B","C","D","E","F","G")))) ))

file_names <- as.list(paste0("T26_05292021/", c("1-NT.CSV", "2-NT.CSV", "1-T.CSV", "2-T.CSV")))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

#plot_2microplates(microplates[[2]], microplates[[4]])

system("more T26_05292021/readme.info") # plate map

s <- data.table(plate=rep(1:4,each=48),
                id=rep(c(15:22, 27:34, 39:46, 51:58, 63:70, 75:82), 4) )

strainDT <- as.data.table(matrix(c(
"classII-3", "YJM2+4-2",  "CLI1-2",    "CLI2-3",    "IFT2-2",    "YCR007C-1", "PRM9-3",    "MST28-1",
"CLI1-4",    "YAR028W-1", "PRM9-4",    "YHL044W-1", "YJM1-3",    "CLI4-1",    "MST28-2",   "YAR027W-3",
"YJM1-1",    "YAR029W-1", "PWF1-1",    "classII-4", "YJM1-2",    "YAR028W-4", "YAR023C-1", "YJM4-4",
"IFT2-1",    "YJM1-4",    "CLI2-1",    "PRM9-1",    "IFT1-2",    "IFT1-1",    "CLI2-4",    "IFT2-4",
"YAR029W-3", "CLI3-4",    "CLI4-3",    "YAR027W-1", "CLI3-3",    "CLI4-4",    "YJM2-4",    "YJM4-3",
"YHL044W-4", "PWF1-3",    "YPS1-2",    "YAR023C-3", "IFT1-3",    "YAR029W-4", "YJM3-3",    "CLI5-2",
"YAR023C-4", "YPS1-1",    "CLI1-3",    "CLI5-3",    "YJM2+4-1",  "CLI4-2",    "IFT2-3",    "YAR029W-2",
"YAR027W-2", "YJM5-2",    "YCR007C-2", "IFT1-4",    "MST28-3",   "YJM2+4-3",  "CLI5-4",    "YJM3-1",
"YJM5-1",    "CLI2-2",    "YJM2-1",    "YJM3-2",    "YJM4-1",    "YAR028W-3", "PRM9-2",    "classII-2",
"YJM4-2",    "YJM2-2",    "YJM2+4-4",  "YHL044W-3", "PWF1-4",    "YAR028W-2", "classII-1", "YJM5-4",
"CLI3-1",    "YCR007C-4", "CLI1-1",    "YPS1-3",    "PWF1-2",    "YJM3-4",    "YAR023C-2", "MST28-4",  
"YPS1-4",    "YJM2-3",    "YCR007C-3", "CLI5-1",    "YAR027W-4", "YHL044W-2", "CLI3-2",    "YJM5-3"  
),
                                 nrow = 8, byrow = FALSE))

strainDT <- melt(strainDT, measure.vars = 1:ncol(strainDT))
s[, "strain" := rep(strainDT$value, 2)]
setkey(s, strain, plate)

largeDT <- s[rep(seq_len(nrow(s)), each = 350), ]
largeDT[, "time" := rep(1:350, 192) * 0.2102857]
largeDT[, "plate_type" := ifelse(plate <= 2, "NT", "T")]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)]
# without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail

largeDT[, "abs" := abs - min(abs)] # absolute blank correction
setcolorder(largeDT, c("plate", "wid", "plate_type", "strain", "rep", "time", "abs"))

ggplot(largeDT, aes(x = time, y = abs, color = plate_type, shape = rep)) +
  xlab("Time (h)") +
  ylab(bquote("Absorbance (OD600)")) +
  geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
  geom_smooth(span = 0.25, size = 0.35, se = FALSE) +
  scale_color_manual(name = "Condition", values = c("black", "deepskyblue1"), labels = c("-K28", "+K28")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(~strain, scales = "free") +
  scale_x_continuous(breaks = seq(0, 72, 72)) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 1)) +
  coord_cartesian(ylim = c(0, largeDT[, max(abs)])) +
  theme_classic() +
  guides(shape = FALSE, color = FALSE) +
  theme(
    legend.position = "top",
    axis.line = element_line(),
    axis.text = element_blank(),
    strip.background = element_blank()
  )
ggsave("T26_smooth_notext.pdf", width = 6, height = 5.25) # export

# Statistical Analysis
s_stat <- copy(s)
s_stat <- s_stat[!(plate >= 3)]
s_stat[, c("strain", "rep") := tstrsplit(strain, "-")]
pheno_fAUC <- largeDT[, AUC(time, abs), by = .(plate, strain, rep)]
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x) {
  return(x[2] / x[1])
}), by = .(strain, rep), .SDcols = "V1"]
s_stat[, "pheno_fAUC" := pheno_fAUC$V1]

# ANOVA
model <- aov(pheno_fAUC ~ strain, data = s_stat)
summary(model)

# Tukey's post-hoc HSD test
TukeyHSD(model)

# Dunnett's post-hoc test with sensitive YAR027W as negative control
DunnettTest(x = s_stat[, pheno_fAUC], g = as.factor(s_stat[, strain]), control = "YAR027W")
