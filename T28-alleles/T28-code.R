# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script.dir <- getwd()
source("../utils.R")

#MAIN SCRIPT

#12/09/2020 T28 - YAR028W Alleles
set.seed(873563)
vec <- sample(paste0(c(rep(c("25","26","27","29","32","33","34","36","37","38"), 3)), c(rep("-1",10), rep("-2",10), rep("-3",10), rep("-4",10))), 40, replace=FALSE)
t(matrix(vec, nrow=10, dimnames = list(c(paste0("Column ", 2:11)), c(paste0("Row ", c("C", "D", "E", "F"))) )))

#       Column 2 Column 3 Column 4 Column 5 Column 6 Column 7 Column 8 Column 9 Column 10 Column 11
# Row C "38-2"   "34-2"   "38-1"   "36-2"   "27-1"   "34-3"   "37-1"   "33-2"   "25-4"    "37-2"   
# Row D "38-3"   "33-3"   "34-4"   "25-2"   "29-4"   "26-3"   "37-3"   "26-1"   "26-2"    "38-4"   
# Row E "32-1"   "36-1"   "32-4"   "29-1"   "27-2"   "32-3"   "36-4"   "27-3"   "33-4"    "29-3"   
# Row F "29-2"   "36-3"   "25-3"   "27-4"   "37-4"   "33-1"   "26-4"   "34-1"   "32-2"    "25-1"  

# blank <- mean(c(.685,.676,.662,.631))
# od <- t(matrix(10*(c(.780,.775,.773,.743,.781,.787,.731,.733,.809,.819,.685,
#                  .625,.730,.795,.755,.715,.804,.695,.782,.747,.829,.676,
#                  .736,.768,.779,.743,.798,.786,.680,.715,.740,.846,.662,
#                  .816,.730,.757,.791,.825,.829,.768,.808,.794,.741,.631) - blank), nrow = 11)) 
# OD ignored - used 2.5 uL inoculum instead


file_names <- as.list(paste0("T28_12092020/", c("NT.CSV", "5050.CSV", "T.CSV")))

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}


print(microplates[[3]])

system("more T28_12092020/readme.info") # plate map

# Plotting
col25 <- colorRampPalette(c("seagreen","white"))(7)
col26 <- colorRampPalette(c("red","white"))(7)
col27 <- colorRampPalette(c("purple","white"))(7)
col29 <- colorRampPalette(c("orange","white"))(7)
col32 <- colorRampPalette(c("orchid1","white"))(7)
col33 <- colorRampPalette(c("brown","white"))(7)
col34 <- colorRampPalette(c("blue","white"))(7)
col36 <- colorRampPalette(c("orange","white"))(7)
col37 <- colorRampPalette(c("skyblue","white"))(7)
col38 <- colorRampPalette(c("turquoise2","white"))(7)

s <- data.table(plate=rep(3,40),
                id=c(26:35, 38:47, 50:59, 62:71)) 

strainDT <- as.data.table(matrix(c("38-2", "34-2", "38-1", "36-2", "27-1", "34-3", "37-1", "33-2", "25-4", "37-2", 
"38-3", "33-3", "34-4", "25-2", "29-4", "26-3", "37-3", "26-1", "26-2", "38-4",  
"32-1", "36-1", "32-4", "29-1", "27-2", "32-3", "36-4", "27-3", "33-4", "29-3",   
"29-2", "36-3", "25-3", "27-4", "37-4", "33-1", "26-4", "34-1", "32-2", "25-1"), nrow = 10, byrow = FALSE) )
strainDT <- melt(strainDT, measure.vars = 1:ncol(strainDT))
s[, "strain" := strainDT$value]
setkey(s, strain)
s[, "col" := c(col25[1:4], col26[1:4], col27[1:4], col29[1:4], col32[1:4], col33[1:4], col34[1:4], col36[1:4], col37[1:4], col38[1:4]) ]
s[, "ltype" := c( rep("solid",4), rep("dotted",4), rep("dotted",4), rep("solid",4), rep("dotted",4), rep("solid",4), rep("dashed",4), rep("dashed",4), rep("dashed",4), rep("solid",4) )]

# Plot growth curves
pattern_ <- "25|38"
curves <- plot_growth_curves(microplates[[1]]$data, ids=c(12,24), col=c("red","green"), title=paste0(pattern_," (T)"))
whichToPlot <- which(s[, strain] %like% pattern_)
min_abs <- blank_correct(s[whichToPlot])
for (p in whichToPlot) {
  curves <- curves + geom_line(mapping = aes_string(x = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$time*0.153394, y = microplates[[ s[p,plate] ]]$data[id == s[p,id] ]$abs - min_abs), color = s[p,col], linetype = s[p,ltype], alpha = 0.7 )
} #.153394 factor = conversion between cycles and hours
( curves <- curves + coord_cartesian(xlim=c(0,42), ylim=c(0,1.2)) ) # plot
#ggsave(paste0("~/Desktop/T30_MSY", pattern_, ".png"), width = 6, height = 4) # export


# Statistical analysis
vec_ttest <- vector("numeric", 0)
for (p in whichToPlot) {
  vec_ttest <- append(vec_ttest, coef(microplates[[ s[p,plate] ]]$nls_list[[ s[p,id] ]])["t0"] )
}
t.test(vec_ttest[1:4], vec_ttest[5:8])


