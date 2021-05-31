library(DescTools)
library(RColorBrewer)

# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd() # is used in code - need to Source script, not Run
source("../utils.R")


#MAIN SCRIPT
# Plate layout
set.seed(85995)
vec <- sample(c(rep("B", 8), paste_expand("18.2-", 1:4), paste_expand(0:20, "-", 1:4)), 96, replace=FALSE)
t(matrix(vec, nrow=8, dimnames = list(c(paste0("Column ", 3:10)), c(paste_expand("Plate ", 1:2, " ", c("B","C","D","E","F","G")))) ))

#          Column 3 Column 4 Column 5 Column 6 Column 7 Column 8 Column 9 Column 10
# Plate 1 B  "17-2"   "8-4"    "20-3"   "18-1"   "15-1"   "3-3"    "19-1"   "10-3"   
# Plate 1 C  "11-3"   "6-4"    "19-4"   "16-3"   "B"      "11-4"   "B"      "4-2"    
# Plate 1 D  "1-2"    "18.2-1" "18.2-4" "1-3"    "0-3"    "13-4"   "14-4"   "2-1"    
# Plate 1 E  "17-4"   "4-1"    "8-3"    "20-4"   "4-3"    "2-3"    "9-1"    "15-4"   
# Plate 1 F  "10-1"   "5-4"    "B"      "B"      "7-4"    "13-1"   "12-4"   "9-4"    
# Plate 1 G  "14-1"   "1-1"    "2-4"    "B"      "B"      "6-1"    "18-3"   "18-2"   

# Plate 2 B  "8-1"    "5-2"    "17-1"   "9-2"    "11-2"   "B"      "3-4"    "7-3"    
# Plate 2 C  "10-2"   "9-3"    "2-2"    "16-4"   "18.2-2" "19-2"   "5-1"    "20-2"   
# Plate 2 D  "12-1"   "7-1"    "6-2"    "16-1"   "18.2-3" "B"      "4-4"    "3-2"    
# Plate 2 E  "0-4"    "13-2"   "12-2"   "18-4"   "12-3"   "8-2"    "19-3"   "20-1"   
# Plate 2 F  "16-2"   "15-2"   "14-2"   "15-3"   "1-4"    "3-1"    "10-4"   "7-2"    
# Plate 2 G  "11-1"   "5-3"    "17-3"   "0-1"    "0-2"    "13-3"   "14-3"   "6-3"    


# T30 - Chimeragenesis (05/06/2021)
file_names <- as.list(paste0("T30_05062021/", c("1-NT.CSV", "2-NT.CSV", "1-T.CSV", "2-T.CSV")))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

system("more T30_05062021/readme.info") # plate map

# # Plotting
s <- data.table(plate=rep(1:4,each=48),
                id=rep(c(15:22, 27:34, 39:46, 51:58, 63:70, 75:82), 4) )

strainDT <- as.data.table(matrix(c("17-2",   "8-4",    "20-3",   "18-1",   "15-1",   "3-3",    "19-1",   "10-3",   
                                   "11-3",   "6-4",    "19-4",   "16-3",   "B",      "11-4",   "B",      "4-2",    
                                   "1-2",    "18.2-1", "18.2-4", "1-3",    "0-3",    "13-4",   "14-4",   "2-1",    
                                   "17-4",   "4-1",    "8-3",    "20-4",   "4-3",    "2-3",    "9-1",    "15-4",   
                                   "10-1",   "5-4",    "B",      "B",      "7-4",    "13-1",   "12-4",   "9-4",    
                                   "14-1",   "1-1",    "2-4",    "B",      "B",      "6-1",    "18-3",   "18-2",   
                                   "8-1",    "5-2",    "17-1",   "9-2",    "11-2",   "B",      "3-4",    "7-3",    
                                   "10-2",   "9-3",    "2-2",    "16-4",   "18.2-2", "19-2",   "5-1",    "20-2",   
                                   "12-1",   "7-1",    "6-2",    "16-1",   "18.2-3", "B",      "4-4",    "3-2",    
                                   "0-4",    "13-2",   "12-2",   "18-4",   "12-3",   "8-2",    "19-3",   "20-1",   
                                   "16-2",   "15-2",   "14-2",   "15-3",   "1-4",    "3-1",    "10-4",   "7-2",    
                                   "11-1",   "5-3",    "17-3",   "0-1",    "0-2",    "13-3",   "14-3",   "6-3"    ),
                                 nrow = 8, byrow = FALSE))

strainDT <- melt(strainDT, measure.vars = 1:ncol(strainDT))
s[, "strain" := rep(strainDT$value, 2)]
setkey(s, strain, plate)

s[which(s[, strain] %in% paste_expand("1-",1:4))]

s[, "col" := rep("red", 96*2)]
s[which(s[, strain] %in% paste_expand("1-",1:4)), "col" := rep("royalblue3", 8)]
s[, "ltype" := rep(c("solid", "dashed"), 96)]

blank_correction <- 0.6954125 # = largeDT[strain == "B", .SD[1:20], wid][, mean(abs)] #check the largeDT code block below

# #Plots
# for (strain in c(0:20, "18.2", "B")) {
#   pattern_list <- c(paste_expand("1-",1:4), paste_expand(strain,"-",1:4))
#   curves <- plot_growth_curves(microplates[[1]]$data, ids=c(12,24), col=c("red","green"), title=strain) #blank plot
#   whichToPlot <- which(s[, strain] %in% pattern_list) # used to be %like% pattern_
#   #min_abs <- blank_correct(s[whichToPlot])
#   for (p in whichToPlot) {
#     curves <- curves + geom_line(mapping = aes_string(x = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$time*0.21104762, y = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$abs - blank_correction), color = s[p,col], linetype = s[p,ltype], alpha = 0.7 )
#   } #0.21104762 factor = conversion between cycles and hours
#   curves <- curves + coord_cartesian(xlim=c(0,67), ylim=c(0,1.1))
#   ggsave(paste0("~/Desktop/T30/chimera_", strain, ".png"), width = 3.5, height = 2) # export
# }


largeDT <- s[rep(seq_len(nrow(s)), each=350), ]
largeDT[, 'time' := rep(1:350, 192)*0.21104762]
largeDT[, 'plate_type' := ifelse(plate <= 2, 'NT', 'T')]
largeDT[, c('strain', 'rep') := tstrsplit(strain, '-')]; names(largeDT)[2] <- 'wid'
largeDT[, 'abs' := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)] #without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, 'abs' := abs - largeDT[strain == "B", .SD[1:20], wid][, mean(abs)] ] # blank correction using average of first 20 measurements across all 8 blank wells
setcolorder(largeDT, c('plate', 'wid', 'plate_type', 'strain', 'rep', 'time', 'abs', 'col'))


ggplot(largeDT[(strain != "B"),], aes(x=time, y=abs, color=plate_type, shape=rep)) +
  xlab("Time (h)") + ylab(bquote('Absorbance (OD600)')) +
  geom_point(mapping=aes(x=time, y=abs), size=0.35, alpha=0 ) +
  geom_smooth(span=0.25, size=0.35, se=FALSE) +
  #geom_line(size=0.25) +
  scale_color_manual(name="Condition", values=c("black", "deepskyblue1"), labels=c("-K28","+K28")) +
  guides(colour=guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(vars(strain), scales = 'free') +
  scale_x_continuous(breaks=seq(0,72,72)) +
  scale_y_continuous(breaks=seq(0.0,1.0,1)) +
  coord_cartesian( ylim=c(0, largeDT[, max(abs)]) ) +
  theme_classic() +
  guides(shape = FALSE, color = FALSE) +
  theme(#axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    #panel.grid.major.x = element_line(),
    legend.position="top",
    axis.line=element_line(),
    axis.text=element_blank(),
    strip.background = element_blank()
  )

ggsave("~/Desktop/T30_smooth_notext.pdf", width = 8, height = 6) # export

s_stat <- copy(s); s_stat <- s_stat[!(plate >= 3 | strain %like% "B")]; s_stat[, 'col' := NULL]
s_stat[, c('strain', 'rep') := tstrsplit(strain, '-')]

pheno_fAUC <- largeDT[!(strain %like% "B"), AUC(time, abs), by = .(plate, strain, rep)]
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x){return( x[2]/x[1] )}), by=.(strain, rep), .SDcols = "V1"]
s_stat[, 'pheno_fAUC' := pheno_fAUC[,3] ]

model <- aov(pheno_fAUC ~ strain, data = s_stat)
summary(model)
TukeyHSD(model)

DunnettTest(x = s_stat[, pheno_fAUC], g = as.factor(s_stat[, strain]), control = "1")



