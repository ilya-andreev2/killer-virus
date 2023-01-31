library(DescTools)
library(RColorBrewer)
library(ggpubr)

# Sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd() # is used in code - need to Source script, not Run
source("../utils.R")


# MAIN SCRIPT
# Plate layout
set.seed(59895)
vec <- sample(c(rep('B', 16), paste_expand(c('0', '1', '2', '3.2', '4', '5', '6', '7', '8', '9.2', '10', '11', '12.2', '13', '14', '15.2', '15.3', '16', '17', '18.2'), '-', 1:4)), 96, replace = FALSE)
t(matrix(c(vec), nrow = 8, dimnames = list(c(paste0("Column ", 3:10)), c(paste_expand("Plate ", 1:2, " ", c("B", "C", "D", "E", "F", "G"))))))

#             Column 3 Column 4 Column 5 Column 6 Column 7 Column 8 Column 9 Column 10
# Plate 1 B   "0-2"    "13-4"   "13-3"   "B"      "8-4"    "16-4"   "17-1"   "7-1"    
# Plate 1 C   "3.2-4"  "15.2-2" "B"      "15.3-3" "4-4"    "B"      "1-3"    "B"      
# Plate 1 D   "3.2-2"  "12.2-2" "16-2"   "5-3"    "B"      "10-3"   "5-1"    "15.3-2" 
# Plate 1 E   "15.2-1" "17-2"   "1-1"    "2-3"    "15.2-4" "7-2"    "15.2-3" "B"      
# Plate 1 F   "8-3"    "12.2-3" "16-1"   "11-4"   "13-2"   "9.2-2"  "7-3"    "0-4"    
# Plate 1 G   "18.2-1" "18.2-4" "17-3"   "10-4"   "B"      "6-2"    "2-1"    "B"      

# Plate 2 B   "7-4"    "B"      "10-1"   "5-2"    "4-1"    "6-4"    "4-3"    "B"      
# Plate 2 C   "3.2-1"  "B"      "11-3"   "1-2"    "B"      "5-4"    "14-1"   "1-4"    
# Plate 2 D   "B"      "14-2"   "4-2"    "13-1"   "11-2"   "12.2-4" "8-2"    "0-3"    
# Plate 2 E   "11-1"   "2-4"    "17-4"   "0-1"    "3.2-3"  "9.2-4"  "9.2-1"  "9.2-3"  
# Plate 2 F   "10-2"   "2-2"    "14-3"   "12.2-1" "B"      "15.3-1" "18.2-2" "B"      
# Plate 2 G   "16-3"   "18.2-3" "14-4"   "6-1"    "15.3-4" "B"      "8-1"    "6-3" 


# T30 - Chimeragenesis (08/09/2021)
file_names <- as.list(paste0("T30_08092021/", c("1-NT.CSV", "2-NT.CSV", "1-T.CSV", "2-T.CSV")))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

s <- data.table(
  plate = rep(1:4, each = 48),
  id = rep(c(15:22, 27:34, 39:46, 51:58, 63:70, 75:82), 4)
)
strainDT <- as.data.table(matrix(
  c(
"0-2",    "13-4",   "13-3",   "B",      "8-4",    "16-4",   "17-1",   "7-1",    
"3.2-4",  "15.2-2", "B",      "15.3-3", "4-4",    "B",      "1-3",    "B",      
"3.2-2",  "12.2-2", "16-2",   "5-3",    "B",      "10-3",   "5-1",    "15.3-2", 
"15.2-1", "17-2",   "1-1",    "2-3",    "15.2-4", "7-2",    "15.2-3", "B",      
"8-3",    "12.2-3", "16-1",   "11-4",   "13-2",   "9.2-2",  "7-3",    "0-4",    
"18.2-1", "18.2-4", "17-3",   "10-4",   "B",      "6-2",    "2-1",    "B",      
"7-4",    "B",      "10-1",   "5-2",    "4-1",    "6-4",    "4-3",    "B",      
"3.2-1",  "B",      "11-3",   "1-2",    "B",      "5-4",    "14-1",   "1-4",    
"B",      "14-2",   "4-2",    "13-1",   "11-2",   "12.2-4", "8-2",    "0-3",    
"11-1",   "2-4",    "17-4",   "0-1",    "3.2-3",  "9.2-4",  "9.2-1",  "9.2-3",  
"10-2",   "2-2",    "14-3",   "12.2-1", "B",      "15.3-1", "18.2-2", "B",      
"16-3",   "18.2-3", "14-4",   "6-1",    "15.3-4", "B",      "8-1",    "6-3" 
    ), nrow = 8, byrow = FALSE))
strainDT <- melt(strainDT, measure.vars = 1:ncol(strainDT))
s[, "strain" := rep(strainDT$value, 2)]
setkey(s, strain, plate)
s[which(s[, strain] %in% paste_expand("1-", 1:4))]
s[, "col" := rep("red", 96 * 2)]
s[which(s[, strain] %in% paste_expand("1-", 1:4)), "col" := rep("royalblue3", 8)]
s[, "ltype" := rep(c("solid", "dashed"), 96)]

blank_correction <- 0.6954125 # = largeDT[strain == "B", .SD[1:20], wid][, mean(abs)] #check the largeDT code block below

largeDT <- s[rep(seq_len(nrow(s)), each = 375), ]
largeDT[, "time" := rep(1:375, 192) * 0.2047111]
largeDT[, "plate_type" := ifelse(plate <= 2, "NT", "T")]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)] # without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, "abs" := abs - largeDT[strain == "B", .SD[1:20], wid][, mean(abs)]] # blank correction using average of first 20 measurements across all 8 blank wells
largeDT <- largeDT[time < 50] #subset to the first 42 hours of measurements, for consistency across similar experiments
setcolorder(largeDT, c("plate", "wid", "plate_type", "strain", "rep", "time", "abs", "col"))

strain_names <- c(
  `B` = "(blank well)", `0` = "(empty vector)", `1` = "KTD1", `2` = "K/U-1", `3.2` = "K/U-2",
  `4` = "K/U-3", `5` = "K/U-4", `6` = "K/U-5", `7` = "K/U-6", `8` = "K/U-7", `9.2` = "9.2",
  `10` = "UIP3", `11` = "U/K-1", `12.2` = "U/K-2", `13` = "U/K-3", `14` = "U/K-4",
  `15.2` = "U/K-5", `15.3` = "U/K-5a", `16` = "U/K-6", `17` = "U/K-7", `18.2` = "18.2"
)
# strains 9.2, and 18.2 not included in downstream analysis

# Figure S11
ggplot(largeDT[!(strain %in% c("9.2", "18.2")), ], aes(x = time, y = abs, color = plate_type, shape = rep)) +
  xlab("Time (h)") +
  ylab(bquote("Absorbance (OD600)")) +
  geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
  geom_smooth(span = 0.50, size = 0.35, se = FALSE) +
  scale_color_manual(name = "Condition", values = c("black", "deepskyblue1"), labels = c("-K28", "+K28")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(vars(strain), labeller = as_labeller(strain_names), scales = "free") +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 1.0)) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1)) +
  theme_classic() +
  guides(shape = "none", color = "none") +
  theme(
    legend.position = "top",
    axis.line = element_line(),
    axis.text = element_blank(),
    strip.background = element_blank()
  )
ggsave("T30_growthCurves.pdf", width = 7, height = 6)

# Statistical analysis
s_stat <- copy(s)
s_stat <- s_stat[!(plate >= 3 | strain %like% "B")]
s_stat[, "col" := NULL]
s_stat[, c("strain", "rep") := tstrsplit(strain, "-")]

pheno_fAUC <- largeDT[!(strain %like% "B"), AUC(time, abs), by = .(plate, strain, rep)]
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x) {
  return(x[2] / x[1])
}), by = .(strain, rep), .SDcols = "V1"]
s_stat[, "pheno_fAUC" := pheno_fAUC[, 3]]

# exclude strains that did not display consistent/replicable behavior from analysis
s_stat <- s_stat[!(strain %in% c("0", "9.2", "18.2"))]

# ANOVA
model <- aov(pheno_fAUC ~ strain, data = s_stat)
summary(model)

# Tukey's post-hoc HSD test
TukeyHSD(model)

tukeyDT <- as.data.table(TukeyHSD(model)$strain, keep.rownames = TRUE); colnames(tukeyDT)[5] <- "p"
tukeyDT[, c("group2", "group1") := tstrsplit(rn, "-")]
tukeyDT[, "p.signif.code" := ifelse(p < 0.001, "***", "*")]
tukeyDT[, "p.e" := sprintf("%1.2e", p)]
tukeyDT[, "p.f" := sprintf("%1.4f", p)]
colnames(tukeyDT) <- c("comp", "diff", "lwr.ci", "upr.ci", "p.adj", "group2", "group1", "p.signif.code", "p.e", "p.f")
setkey(tukeyDT, p.adj)

# Compare Chimera 2 (K/U-1) to every other strain
tukeyDT[group1 == "2"]
tukeyDT[group2 == "2"]

# Compare everything to Chimera 10 (UIP3)
tukeyDT[group1 == "10"]
tukeyDT[group2 == "10"]

# Figure 5A
ggstripchart(s_stat[!(strain %in% c("0", "9.2", "15.3", "18.2"))], "strain", "pheno_fAUC",
             xlab = "Strain", ylab = "K28 Resistance",
             add = "mean", add.params = list(color = "red"),
             # order = rev(c("0", "1", "2", "3.2", "4", "5", "6", "7", "8", "10", "11", "12.2", "13", "14", "15.2", "16", "17"))
             order = rev(c("10", "11", "12.2", "13", "14", "15.2", "16", "17", 
                           "1", "2", "3.2", "4", "5", "6", "7", "8"))
             ) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.5)) +
  coord_flip()
ggsave("T30_tall.pdf", width = 4, height = 11)

