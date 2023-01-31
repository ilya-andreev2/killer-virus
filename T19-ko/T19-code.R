# Sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")

#MAIN SCRIPT
file_names <- paste0("T19_02282020/", list("TRno152.CSV"))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

# Plate Layout
# Row A: MSY1 (BY MATa)
# Row B: KO-YAR027 (BY MATa)
# Row C: KO-YAR028 (BY MATa)
# Row D: MSY52 (diploid M28-infected)
# Rows E,F: MSY31
# Rows G,H: MSY32

# Percent toxin by column: 
seq(160, 0, length.out = 10)/160*100

s <- data.table( plate=rep(rep(1, each=8), 2),
                 id=c(12*0:7+2, 12*0:7+3, 12*0:7+4, 12*0:7+5, 12*0:7+6, 12*0:7+7,
                      12*0:7+8, 12*0:7+9, 12*0:7+10, 12*0:7+11) #2nd column ids T -> 11th for NT
)

strainDT <- as.data.table(matrix(c(
  "MSY1-1", "YAR027Wdel-1", "YAR028Wdel-1", "MSY52-1", "MSY31-1", "MSY31-2", "MSY32-1", "MSY32-2"), 
  nrow = 8, byrow = TRUE) )
strainDT <- melt(strainDT, measure.vars = 1:ncol(strainDT))
s[, "strain" := rep(strainDT$value, 10)]
setkey(s, strain, id)
s[, "type" := rep(c("T","T2","T3","T4","T5","T6","T7","T8","T9","NT"), 8)]

largeDT <- s[rep(seq_len(nrow(s)), each = 300), ]
largeDT[, "time" := rep(1:300, 80) * 0.1323333]
largeDT[, c("strain", "rep") := tstrsplit(strain, "-")]
names(largeDT)[2] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)]
# largeDT[, "abs" := abs - min(abs)] # absolute blank correction
largeDT[, "abs" := abs - largeDT[(strain == "YAR028Wdel") & (type == "T"), .SD[1:10], wid][, mean(abs)]] # blank correction using average of first 10 measurements across all wells
setcolorder(largeDT, c("plate", "wid", "type", "strain", "rep", "time", "abs"))

# Figure S4
ggplot(largeDT[strain %in% c("MSY1", "YAR027Wdel", "YAR028Wdel")], aes(x = time, y = abs, color = type, shape = rep)) +
  xlab("Time (h)") +
  ylab(bquote("Absorbance (OD600)")) +
  geom_point(mapping = aes(x = time, y = abs), size = 0.35, alpha = 0) +
  geom_smooth(span = 0.35, size = 0.35, se = FALSE) +
  scale_color_manual(name = "Condition", values = c("#000000", colorRampPalette(c("deepskyblue1","black"))(15)[1:9]), labels = c("-K28", "100% K28", "88.9% K28", "77.8% K28", "66.7% K28", "55.5% K28", "44.4% K28", "33.3% K28", "22.2% K28", "11.1% K28")) +
  # conditions get sorted first - hence why the order is weird
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(vars(strain), scales = "free") +
  scale_x_continuous(breaks = c(0, 39)) +
  scale_y_continuous(breaks = seq(0.0, 1.5, 0.5)) +
  coord_cartesian(ylim = c(0, largeDT[, max(abs)])) +
  theme_classic() +
  guides(shape = "none") +
  theme(
    #legend.position = "top",
    axis.line = element_line(),
    axis.text = element_blank(),
    strip.background = element_blank()
  )
ggsave("T19_growthCurves.pdf", width = 8, height = 3)
