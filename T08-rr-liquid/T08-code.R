library(DescTools)
library(ggpubr)

# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")


# MAIN SCRIPT

# 04/09/2021 - T08 Round 3 (4xBiol Rep + 5xADE to MSY39)
# same as 03/18/2021, but stronger 48h toxic supernatant
# plate map and additional info is contained in "T08_04092021/readme.info" file
set.seed(31112)
vec <- sample(c(paste_expand("STD-", 1:16), paste_expand(24:39, "-", 1:5)), 96, replace = FALSE) #randomized vector of sample names. 
# Format: "XX-Y", where XX refers to MSY XX strain number, Y denotes the biological replicate index, and STD wells are for generating the Toxin Standard Curve

t(matrix(vec,
  nrow = 8,
  dimnames = list(
    c(paste0("Column ", 3:10)),
    c(paste_expand("Plate ", 1:2, " ", c("B", "C", "D", "E", "F", "G"))
    )
  )
)) # assign well positions

# compute Toxin:Non-toxin ratio of supernatants for STD-1 through STD-16 (full toxic to full non-toxic, respectively)
vec_std <- paste0( 
  sprintf("%.1f", round(seq(160, 0, length.out = 16), 1)), ":",
  sprintf("%.1f", round(seq(0, 160, length.out = 16), 1), 1)
)

# Convert above to a named matrix
mat_std <- matrix(vec_std,
  nrow = 1, byrow = TRUE,
  dimnames = list(NULL, paste_expand("STD-", 1:16))
)

file_names <- as.list(paste0("T08_04092021/", c("1-NT.CSV", "2-NT.CSV", "1-T.CSV", "2-T.CSV")))

microplates <- list()
for (i in seq_len(length(file_names))) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

# Example plot of the first microplate
microplates[[1]]

# Processing
s <- data.table(
  plate = rep(1:4, each = 48),
  id = rep(c(15:22, 27:34, 39:46, 51:58, 63:70, 75:82), 4)
)

strainDT <- as.data.table(matrix(c(
  "33-4", "31-5", "STD-6", "37-2", "32-2", "30-3", "29-2", "26-5",
  "25-1", "35-3", "32-4", "24-1", "STD-7", "STD-13", "37-1", "35-2",
  "37-4", "26-4", "31-2", "24-5", "STD-5", "34-2", "32-1", "33-2",
  "STD-15", "31-3", "STD-12", "STD-10", "39-2", "36-5", "27-1", "24-4",
  "30-1", "30-2", "34-1", "28-4", "38-4", "35-1", "STD-3", "36-4",
  "24-3", "28-1", "24-2", "25-3", "31-1", "33-5", "28-3", "39-3",
  "27-3", "34-3", "25-5", "38-3", "26-2", "37-3", "STD-2", "STD-9",
  "27-2", "STD-14", "29-1", "36-2", "39-4", "29-3", "36-1", "28-5",
  "28-2", "33-1", "25-4", "27-4", "30-5", "26-3", "STD-1", "25-2",
  "29-5", "38-1", "30-4", "35-4", "34-5", "39-5", "38-2", "STD-4",
  "35-5", "34-4", "36-3", "32-3", "37-5", "STD-11", "29-4", "STD-16",
  "33-3", "32-5", "STD-8", "39-1", "38-5", "27-5", "31-4", "26-1"
),
nrow = 8, byrow = FALSE
))

strainDT <- melt(strainDT, measure.vars = seq_len(ncol(strainDT)))
s[, "strain" := rep(strainDT$value, 2)]
setkey(s, strain, plate)
s[, "col" := ifelse(plate <= 2, "black", "royalblue3")]
s[, "plot_index" := match(
  tstrsplit(strain, "-")[[1]],
  c(
    "26", "24", "37", "33", 
    "29", "36", "35", "25",
    "31", "28", "34", "39", 
    "32", "38", "27", "30"
  )
)]

largeDT <- s[rep(seq_len(nrow(s)), each = 375), ]
largeDT[, "time" := rep(1:375, 192) * 0.2085333]
largeDT[, "plate_type" := ifelse(plate <= 2, "NT", "T")]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)]
# without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, "abs" := abs - min(abs)] # absolute blank correction
largeDT <- largeDT[time < 42] #subset to the first 42 hours of measurements, for consistency across different experiments
setcolorder(largeDT, c("plate", "wid", "plate_type", "strain", "rep", "time", "abs", "col"))

strain_names <- c(
  `1` = "RM", `2` = "M22", `3` = "YJM981", `4` = "I14",
  `5` = "CLIB413", `6` = "273614", `7` = "PW5", `8` = "BY",
  `9` = "YJM454", `10` = "YJM145", `11` = "Y10", `12` = "CLIB219",
  `13` = "YPS1009", `14` = "CBS2888", `15` = "YPS163", `16` = "YJM978"
)

# Figure 1B
ggplot(largeDT[(strain != "STD"), ], aes(x = time, y = abs, color = plate_type, shape = rep)) +
  xlab("Time (h)") +
  ylab(bquote("Absorbance (OD600)")) +
  geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
  geom_smooth(span = 0.25, size = 0.35, se = FALSE) +
  scale_color_manual(name = "Condition", values = c("black", "deepskyblue1"), labels = c("-K28", "+K28")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(vars(plot_index), labeller = as_labeller(strain_names), scales = "free") +
  scale_x_continuous(breaks = c(0, 42)) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 1)) +
  coord_cartesian(ylim = c(0, largeDT[, max(abs)])) +
  theme_classic() +
  guides(shape = "none", color = "none") +
  theme(
    legend.position = "top",
    axis.line = element_line(),
    axis.text = element_blank(),
    strip.background = element_blank()
  )
ggsave("T08_growthCurves.pdf", width = 6, height = 5.25) # export


# Statistical analysis
s_stat <- copy(s)
s_stat <- s_stat[!(plate >= 3 | strain %like% "STD")]
s_stat[, "col" := NULL]
s_stat[, c("strain", "rep") := tstrsplit(strain, "-")]
pheno_fAUC <- largeDT[!(strain %like% "STD"), AUC(time, abs), by = .(plate, strain, rep)]
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x) {
  return(x[2] / x[1])
}), by = .(strain, rep), .SDcols = "V1"]
s_stat[, "pheno_fAUC" := pheno_fAUC$V1]

# ANOVA
model <- aov(pheno_fAUC ~ strain, data = s_stat)
summary(model)

# Tukey's post-hoc HSD test
tukeyDT <- as.data.table(TukeyHSD(model)$strain, keep.rownames = TRUE); colnames(tukeyDT)[5] <- "p"
tukeyDT[, c("group2", "group1") := tstrsplit(rn, "-")]
tukeyDT[, "p.signif.code" := ifelse(p < 0.001, "***", "*")] # this really sucks but OK for now
tukeyDT[, "p.e" := sprintf("%1.2e", p)]
colnames(tukeyDT) <- c("comp", "diff", "lwr.ci", "upr.ci", "p.adj", "group2", "group1", "p.signif.code", "p.e")

# Dunnett's post-hoc test with resistant BY as control - not performed
# DunnettTest(x = s_stat[, pheno_fAUC], g = as.factor(s_stat[, strain]), control = "25")

# Figure S1
ggstripchart(s_stat, "strain", "pheno_fAUC",
             xlab = "Strain", ylab = "K28 Resistance",
             add = "mean", add.params = list(color = "red"),
             order =   c(
               "26", "24", "37", "33", 
               "29", "36", "35", "25",
               "31", "28", "34", "39", 
               "32", "38", "27", "30"
             )
) + 
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.5)) +
  stat_pvalue_manual(tukeyDT[comp %like% "25" & p.adj < 0.05], label = "p.e", y.position = seq(1.3, 1.8, 0.12))
ggsave("T08_stats.pdf", width = 8, height = 6) # export
