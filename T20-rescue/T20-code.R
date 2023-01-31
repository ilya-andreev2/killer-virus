library(DescTools)
library(ggpubr)

# Sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")

#MAIN SCRIPT
# 09/04/2020 T20 - YAR028W Phenotypic Rescue Experiment (Fig. 2 D and E)
set.seed(1029)
vec <- sample(paste0("MSY", rep(as.character(91:99), 3), c(rep("-1",9), rep("-2",9), rep("-3",9), rep("-4",9))), 36, replace=FALSE)
t(matrix(vec, ncol = 6, dimnames = list(c(paste0("Plate ", 1:6)), c(paste0("Row ", c("B","C","D","E","F","G"))))))

# this experiment is different, in that each row on an individual plate contains the same
# strain, which is spread across all columns. T supernatant in 10 gradations, 
# from left (all Toxic) to right (all Non-Toxic)

file_names <- paste0("T20_09042020/", list("Plate1.CSV","Plate2.CSV","Plate3.CSV","Plate4.CSV","Plate5.CSV","Plate6.CSV","Plate7.CSV") )

microplates <- list()
for (i in 1:length(file_names)) {
   microplates[[i]] <- process_microplate(file_names[[i]])
}

# Plotting
colv91 <- colorRampPalette(c("seagreen","white"))(7)
colv92 <- colorRampPalette(c("red","white"))(7)
colv93 <- colorRampPalette(c("orange","white"))(7)
colv94 <- colorRampPalette(c("#FFC000","white"))(7)
colv95 <- colorRampPalette(c("orchid1","white"))(7)
colv96 <- colorRampPalette(c("#999933","white"))(7)
colv97 <- colorRampPalette(c("#0070C0","white"))(7)
colv98 <- colorRampPalette(c("#FFC000","white"))(7)
colv99 <- colorRampPalette(c("#999933","white"))(7)

# TOXIC (+2) and NON-TOXIC (+11) curves
s <- data.table( plate=rep(rep(1:6, each=6), 2),
                      id=c(rep(12*1:6+2,6), rep(12*1:6+11,6)) #2nd column ids -> 11th for NT
)

strainDT <- as.data.table(matrix(c(
          "MSY91-2", "MSY98-4", "MSY96-1", "MSY95-2", "MSY98-1", "MSY97-3",
          "MSY93-3", "MSY96-3", "MSY94-2", "MSY95-3", "MSY92-3", "MSY97-2",
          "MSY99-3", "MSY91-4", "MSY94-1", "MSY99-2", "MSY98-3", "MSY92-2",
          "MSY94-4", "MSY93-2", "MSY96-2", "MSY97-1", "MSY95-1", "MSY97-4",
          "MSY98-2", "MSY93-1", "MSY91-3", "MSY93-4", "MSY92-4", "MSY99-4",
          "MSY95-4", "MSY92-1", "MSY96-4", "MSY91-1", "MSY94-3", "MSY99-1"), 
          nrow = 6, byrow = TRUE) )
strainDT <- melt(strainDT, measure.vars = 1:ncol(strainDT))
s[, "strain" := rep(strainDT$value, 2)]
setkey(s, strain, id)
s[, "type" := rep(c("T","NT"), 36)]

#needed for non-facet plotting
s[, "col" := rep(c(colv91[1:4], colv92[1:4], colv93[1:4], colv94[1:4], colv95[1:4], colv96[1:4], colv97[1:4], colv98[1:4], colv99[1:4]), each = 2) ]
s[, "ltype" := rep(c( rep("solid",4), rep("dotted",4), rep("dotted",4), rep("solid",4), rep("dotted",4), rep("solid",8), rep("dashed", 8)), each = 2)]

largeDT <- s[rep(seq_len(nrow(s)), each = 187), ]
largeDT[, "time" := rep(1:187, 72) * 0.4277184]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)]
# without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, "abs" := abs - largeDT[, .SD[1:10], wid][, mean(abs)]]
#largeDT[, "abs" := abs - min(abs)] # absolute blank correction - alternative to above
largeDT <- largeDT[time < 42] #subset to the first 42 hours of measurements, for consistency across similar experiments
setcolorder(largeDT, c("plate", "wid", "type", "strain", "rep", "time", "abs", "col", "ltype"))

strain_names <- c(
   `MSY93` = "RM[RM-KTD1]", `MSY94` = "RM[BY-KTD1]", `MSY95` = "ktd1D[RM-KTD1]",
   `MSY96` = "ktd1D[BY-KTD1]", `MSY97` = "BY", `MSY98` = "RM", `MSY99` = "ktd1D"
)

# Whole experiment plot
ggplot(largeDT[!(strain %in% c("MSY91", "MSY92")),], aes(x = time, y = abs, color = type, shape = rep)) +
   xlab("Time (h)") +
   ylab(bquote("Absorbance (OD600)")) +
   geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
   geom_smooth(span = 0.25, size = 0.35, se = FALSE) +
   scale_color_manual(name = "Condition", values = c("black", "deepskyblue1"), labels = c("-K28", "+K28")) +
   guides(color = guide_legend(override.aes = list(alpha = 1))) +
   facet_wrap(vars(strain), labeller = as_labeller(strain_names), scales = "free") +
   scale_x_continuous(breaks = c(0, 42)) +
   scale_y_continuous(breaks = seq(0.0, 1.5, 0.5)) +
   coord_cartesian(ylim = c(0, largeDT[, max(abs)])) +
   theme_classic() +
   guides(shape = "none", color = "none") +
   theme(
      legend.position = "top",
      axis.line = element_line(),
      axis.text = element_blank(),
      strip.background = element_blank()
   )
#ggsave("T20_growthCurves.pdf", width = 5, height = 5) # export


pattern_ <- "MSY94|MSY98" #RM[BY-YAR028W]|RM[empty] # uncomment this for Fig. 2E
# pattern_ <- "MSY96|MSY97|MSY99" #yar028wD[BY-YAR028W]|BY[empty]|yar028wD[empty] # uncomment this for Fig. 2D

gc_plot <- plot_growth_curves(microplates[[1]]$data, ids=c(12,24), col=c("red","green"), title="Growth in +K28")#title=paste0(pattern_," (T)"))
whichToPlot <- which(s[type == "T", strain] %like% pattern_)
min_abs <- blank_correct(s[whichToPlot])
for (p in whichToPlot) {
  gc_plot <- gc_plot + 
   #geom_point(mapping = aes_string(x = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$time*0.427718, y = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$abs - min_abs), color = s[p,col], size = 0.3, alpha = 0.5 )  
   geom_smooth(mapping = aes_string(x = microplates[[ s[type == "T"][p,plate] ]]$data[id == s[type == "T"][p,id] ]$time*0.427718, y = microplates[[ s[type == "T"][p,plate] ]]$data[id == s[type == "T"][p,id] ]$abs - min_abs), color = s[type == "T"][p,col], alpha = 0.3,  span = 0.25, size = 0.35, se = FALSE)
  #geom_line(mapping = aes_string(x = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$time*0.427718, y = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$abs - min_abs), color = s[p,col], linetype = s[p,ltype], alpha = 0.7 )
} #.427718 factor = conversion between cycles and hours

( gc_plot <- gc_plot + 
      #geom_smooth(span = 0.25, size = 0.35, se = FALSE) +
      scale_x_continuous(breaks = c(0, 42)) +
      coord_cartesian(xlim=c(0,42), ylim=c(0,1.5)) ) # plot

# ggsave("94vs98_T.pdf", width = 3.5, height = 2.8) #2E
# ggsave("96vs97vs99_T.pdf", width = 3.5, height = 2.8) #2D

# Statistical analysis
s_stat <- copy(s)
s_stat <- s_stat[!(strain %in% paste_expand("MSY9", 1:2, "-", 1:4))] #remove MSY91 and 92 from analysis
s_stat <- s_stat[type == "T"]
s_stat[, c("col", "ltype", "type") := NULL]
s_stat[, c("strain", "rep") := tstrsplit(strain, "-")]

pheno_fAUC <- largeDT[!(strain %in% c("MSY91", "MSY92")), AUC(time, abs), by = .(type, strain, rep)]
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x) {
   return(x[1] / x[2])
}), by = .(strain, rep), .SDcols = "V1"]
s_stat[, "pheno_fAUC" := pheno_fAUC$V1]


# Specific directional t-tests : 97 vs 99 = BY vs ktd1D
# 96 vs 99 = ktd1D[KTD1] vs ktd1D
# 98 vs 94 = RM vs RM[KTD1]
# 93 vs 98 = RM[RM-KTD1] vs RM (no effect)
# 95 vs 99 = ktd1D[RM-KTD1] vs ktd1D
my_comparisons <- list(c("MSY97", "MSY99"), c("MSY96", "MSY99"), c("MSY94", "MSY98"), c("MSY93", "MSY98"), c("MSY95", "MSY99"))

# Figure S6
ggstripchart(s_stat, "strain", "pheno_fAUC", shape = 1,
             xlab = "Strain", ylab = "K28 Resistance",
             add = "mean", add.params = list(color = "red"), 
             #order = c("MSY93", "MSY94", "MSY95", "MSY96", "MSY97", "MSY98", "MSY99")) + 
             order = c("MSY97", "MSY99", "MSY96", "MSY95", "MSY98", "MSY94", "MSY93")) +
   #stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE) # 2-tailed t test
   stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE, method.args = list(alternative = "greater")) #1-tailed t test
ggsave("T20_stats.pdf", width = 6.5, height = 5.2)
