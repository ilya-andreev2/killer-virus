library(DescTools)
library(ggpubr)

# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")


# MAIN SCRIPT
# Plate layout
set.seed(873563) #these odd seeds come from when the universe is speaking to me. dont @ me!
vec <- sample(paste0(c(rep(c("25", "26", "27", "29", "32", "33", "34", "36", "37", "38"), 3)), c(rep("-1", 10), rep("-2", 10), rep("-3", 10), rep("-4", 10))), 40, replace = FALSE)
t(matrix(vec, nrow = 10, dimnames = list(c(paste0("Column ", 2:11)), c(paste0("Row ", c("C", "D", "E", "F"))))))

# 12/09/2020 T28 - YAR028W Alleles in ktd1-KO
file_names <- as.list(paste0("T28_12092020/", c("NT.CSV", "5050.CSV", "T.CSV")))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

s <- data.table(
  plate = rep(c(1,3), each = 40),
  id = c(26:35, 38:47, 50:59, 62:71)
)

strainDT <- as.data.table(matrix(c(
  "38-2", "34-2", "38-1", "36-2", "27-1", "34-3", "37-1", "33-2", "25-4", "37-2",
  "38-3", "33-3", "34-4", "25-2", "29-4", "26-3", "37-3", "26-1", "26-2", "38-4",
  "32-1", "36-1", "32-4", "29-1", "27-2", "32-3", "36-4", "27-3", "33-4", "29-3",
  "29-2", "36-3", "25-3", "27-4", "37-4", "33-1", "26-4", "34-1", "32-2", "25-1"
), nrow = 10, byrow = FALSE))
strainDT <- melt(strainDT, measure.vars = seq_len(ncol(strainDT)))
s[, "strain" := rep(strainDT$value, 2)]
setkey(s, strain, plate)

largeDT <- s[rep(seq_len(nrow(s)), each = 275), ]
largeDT[, "time" := rep(1:275, 80) * 0.153394]
largeDT[, "plate_type" := ifelse(plate == 1, "NT", "T")]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)] # without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, "abs" := abs - largeDT[, .SD[1:20], wid][, mean(abs)]] # blank correction using average of first 20 measurements across all 8 blank wells
setcolorder(largeDT, c("plate", "wid", "plate_type", "strain", "rep", "time", "abs"))

# Figure 3B (including identical alleles)
ggplot(largeDT, aes(x = time, y = abs, color = plate_type, shape = rep)) +
  xlab("Time (h)") +
  ylab(bquote("Absorbance (OD600)")) +
  geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
  geom_smooth(span = 0.35, size = 0.35, se = FALSE) +
  scale_color_manual(name = "Condition", values = c("black", "deepskyblue1"), labels = c("-K28", "+K28")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(vars(strain), scales = "free") +
  scale_x_continuous(breaks = seq(0, 48, 48)) +
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
ggsave("T28_growthCurves.pdf", width = 8, height = 6) # export

# Statistical analysis
s_stat <- copy(s)
s_stat <- s_stat[!(plate == 3)]
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

tukeyDT <- as.data.table(TukeyHSD(model)$strain, keep.rownames = TRUE); colnames(tukeyDT)[5] <- "p"
tukeyDT[, c("group2", "group1") := tstrsplit(rn, "-")]
tukeyDT[, "p.signif.code" := ifelse(p < 0.001, "***", "*")]
tukeyDT[, "p.e" := sprintf("%1.2e", p)]
colnames(tukeyDT) <- c("comp", "diff", "lwr.ci", "upr.ci", "p.adj", "group2", "group1", "p.signif.code", "p.e")
setkey(tukeyDT, p.adj)

# Figure S7B
ggstripchart(s_stat, "strain", "pheno_fAUC",
             xlab = "Strain", ylab = "K28 Resistance",
             add = "mean", add.params = list(color = "red"),
             order =   c(
               "25", "27", "29", "34", "36", 
               "38", "32", "33", "26", "37"
             )
) + 
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.5)) +
  stat_pvalue_manual(tukeyDT[comp %like% "25" & p.adj < 0.05], label = "p.e", y.position = c(1.42, 1.30, 1.54, 0.70))
ggsave("T28_stats.pdf", width = 6.5, height = 5.2)

