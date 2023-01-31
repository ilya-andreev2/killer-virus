library(ggpubr)
library(DescTools)

# Sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")

#MAIN SCRIPT
#06/18/2021 K28 Toxin Adsorption (attempt #3: 15 min adsorption)
set.seed(8686030)
vec <- sample(c(
  paste_expand("STD-", 1:12),
  paste_expand(c("94", "96", "97", "98", "99", "100", "104", "T", "NT"), "-", 1:4)
), 48, replace = FALSE)

t(matrix(vec, nrow=8, dimnames = list(
  c(paste0("Column ", 3:10)), 
  c(paste_expand(
  "Row ",
  c("B", "C", "D", "E", "F", "G")
  )) )))

# ignore STD wells - did not have sufficient T supernatant
#       Column 3 Column 4 Column 5 Column 6 Column 7 Column 8 Column 9 Column 10
# Row B "STD-10" "99-1"   "NT-2"   "STD-4"  "STD-7"  "104-4"  "STD-2"  "98-1"   
# Row C "94-1"   "STD-9"  "NT-1"   "T-3"    "104-3"  "T-4"    "STD-5"  "94-4"   
# Row D "100-1"  "96-2"   "96-4"   "99-4"   "97-4"   "100-4"  "96-3"   "94-3"   
# Row E "98-4"   "T-2"    "T-1"    "STD-6"  "104-2"  "STD-8"  "NT-3"   "STD-3"  
# Row F "98-3"   "100-2"  "96-1"   "98-2"   "STD-1"  "94-2"   "99-3"   "100-3"  
# Row G "97-3"   "104-1"  "97-1"   "NT-4"   "STD-12" "97-2"   "STD-11" "99-2" 

vec_std <- sample(paste0(sprintf("%.1f", round(seq(800, 0, length.out = 12), 1)), ':', sprintf("%.1f", round(seq(0,800, length.out=12), 1), 1)), 12, replace=FALSE)
names(vec_std) <- paste_expand("STD-", 1:12)

# manual OD600 measurements of the 30 mL cultures
t(matrix( 50/c(3.66, 3.66, 3.71, 3.82, #MSY94
               3.73, 3.84, 3.72, 3.64, #MSY96
               2.88, 3.91, 3.77, 2.48, #MSY97
               3.66, 4.00, 3.24, 3.55, #MSY98
               3.61, 3.70, 3.30, 4.18, #MSY99
               3.53, 3.48, 3.50, 3.83, #MSY100
               3.93, 4.25, 4.08, 3.76 #MSY104
               ), nrow = 4))

file_names <- as.list(paste0("T25_06182021/", c("1.CSV", "2.CSV", "3.CSV")))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

s <- data.table(
  plate = rep(1:3, each = 48),
  id = rep(c(15:22, 27:34, 39:46, 51:58, 63:70, 75:82), 3)
)

strainDT <- as.data.table(matrix(c(
  "STD-10", "99-1",   "NT-2",   "STD-4",  "STD-7",  "104-4",  "STD-2",  "98-1",   
  "94-1",   "STD-9",  "NT-1",   "T-3",    "104-3",  "T-4",    "STD-5",  "94-4",   
  "100-1",  "96-2",   "96-4",   "99-4",   "97-4",   "100-4",  "96-3",   "94-3",   
  "98-4",   "T-2",    "T-1" ,   "STD-6",  "104-2",  "STD-8",  "NT-3",   "STD-3",  
  "98-3",   "100-2",  "96-1",   "98-2",   "STD-1",  "94-2",   "99-3",   "100-3",  
  "97-3",   "104-1",  "97-1",   "NT-4",   "STD-12", "97-2",   "STD-11", "99-2" 
),
nrow = 8, byrow = FALSE
))

strainDT <- melt(strainDT, measure.vars = seq_len(ncol(strainDT)))
s[, "strain" := rep(strainDT$value, 3)]
setkey(s, strain, plate)

s[, "plot_index" := match(
  tstrsplit(strain, "-")[[1]],
  c(
    "94", "96", "97", "98", "99", "100", "104", "T", "NT"
  )
)]

largeDT <- s[rep(seq_len(nrow(s)), each = 350), ]
largeDT[, "time" := rep(1:350, 144) * 0.16704762]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)]
largeDT[, "abs" := abs - largeDT[, .SD[10:30], by = c("wid", "plate")][, mean(abs)]] # blank correction using average of first 20 measurements across all 8 blank wells
largeDT <- largeDT[time < 42] #subset to the first 42 hours of measurements, for consistency across similar experiments
setcolorder(largeDT, c("plate", "wid", "strain", "rep", "time", "abs"))

strain_names <- c(
  `100` = "mnn2D+[empty vector]", `104` = "BY+[KTD1BY]", `94` = "RM+[KTD1BY]", `96` = "ktd1D+[KTD1BY]",
  `97` = "BY+[empty vector]", `98` = "RM+[empty vector]", `99` = "ktd1D+[empty vector]",
  `NT` = "(non-toxic cell-free media)", `T` = "(toxic cell-free media)"
)

# Figure S8
ggplot(largeDT[(strain != "STD" & plate == 3), ], aes(x = time, y = abs, shape = rep, color = plate)) +
  xlab("Time (h)") +
  ylab(bquote("Absorbance (OD600)")) +
  geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
  geom_smooth(span = 0.25, size = 0.35, se = FALSE) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(vars(strain), labeller = as_labeller(strain_names), scales = "free") +
  scale_x_continuous(breaks = c(0, 42)) +
  scale_y_continuous(breaks = seq(0.0, 1, 1)) +
  coord_cartesian(ylim = c(0, largeDT[, max(abs)])) +
  theme_classic() +
  guides(shape = "none", color = "none") +
  theme(
    legend.position = "top",
    axis.line = element_line(),
    axis.text = element_blank(),
    strip.background = element_blank()
  )
ggsave("T25_growthCurves.pdf", width = 6, height = 6) # export

# Statistical Analysis
s_stat <- copy(s)
s_stat <- s_stat[!(strain %like% "STD") & plate == 3]
s_stat[, "plot_index" := NULL]
s_stat[, c("strain", "rep") := tstrsplit(strain, "-")]
pheno_AUC <- largeDT[!(strain %like% "STD") & plate == 3, AUC(time, abs), by = .(plate, strain, rep)]
s_stat[, "pheno_AUC" := pheno_AUC$V1]

# ANOVA
model <- aov(pheno_AUC ~ strain, data = s_stat)
summary(model)

# Tukey's post-hoc HSD test
tukeyDT <- as.data.table(TukeyHSD(model)$strain, keep.rownames = TRUE); colnames(tukeyDT)[5] <- "p"
tukeyDT[, c("group2", "group1") := tstrsplit(rn, "-")]
tukeyDT[, "p.signif.code" := ifelse(p < 0.001, "***", "*")] # this kinda sucks but OK for this purpose
tukeyDT[, "p" := sprintf("%1.3f", p)]
colnames(tukeyDT) <- c("comp", "diff", "lwr.ci", "upr.ci", "p", "group2", "group1", "p.signif.code")
tukeyDT

# Dunnett's post-hoc test, with RM growth curve in non-toxic supernatant as control
dtest <- DunnettTest(x = s_stat[, pheno_AUC], g = as.factor(s_stat[, strain]), control = "100")
dtest <- as.data.table(dtest$`100`, keep.rownames = TRUE)
dtest[, c("group2", "group1") := tstrsplit(rn, "-")]
dtest[, "p.signif.code" := ifelse(pval < 0.001, "***", "*")] 
colnames(dtest) <- c("comp", "diff", "lwr.ci", "upr.ci", "p", "group2", "group1", "p.signif.code")
dtest

# Figure 4A
my_comparisons <- list(c("98", "94"), c("99", "96"), c("97", "104"), c("97", "100"), c("98", "97"))
ggstripchart(s_stat[!(strain %in% c("T", "NT"))], "strain", "pheno_AUC", shape = 1,
             xlab = "Strain", ylab = "Toxin Depletion (AUC, arbitrary units)",
             add = "mean", add.params = list(color = "red"), 
             order = c("98", "94", "97", "104", "99", "96", "100")) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE, method.args = list(alternative = "greater")) #1-tailed t test
ggsave("T25_stats.pdf", width = 5, height = 3.5) # export


