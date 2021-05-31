# sets working directory to the script location - only works in RStudio
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

script.dir <- dirname(sys.frame(1)$ofile) # is used in code - need to Source script, not Run
setwd(script.dir)

library(ggplot2)
library(gdata)
library(stringr)
library(data.table)


#OBJECTS
new_microplate <- function(data = data.table(),
                           nls_list = list(),
                           endpoint = FALSE) {
  # <<<<<<<<<TO IMPLEMENT: name, date, !cycle_time!, protocol, other metadata from .csv>>>>>>>>>>
  stopifnot(is.data.table(data))
  plate <- list(data = data, nls_list = nls_list, endpoint = endpoint)
  attr(plate, "class") <- "microplate"
  plate
}

print.microplate <- function(microplate, col="magenta", ...) {
  plt <- ggplot(data = microplate$data) +
    geom_point(mapping=aes(x=time, y=abs), color=col) +
    geom_line(mapping=aes(x=time, y=fit)) +
    facet_wrap(~well, nrow=8) #+
    #coord_cartesian(xlim=c(0,500))
  print(plt)
  print(microplate$data)
}


#FUNCTIONS
fit_nls_curve <- function(nls.model, time.points, na.flag) {
  if(na.flag == TRUE) { return( as.double(rep(NA, length(time.points))) ) }
  else { return( predict(nls.model, newdata=time.points) ) }
}

blank_correct <- function(s) {
  min_abs <- 999999
  for (i in 1:nrow(s)) {
    print(min(microplates[[ s[i,plate] ]]$data[ id == s[i,id] ]$abs))
    min_temp <- min(microplates[[ s[i,plate] ]]$data[ id == s[i,id] ]$abs, na.rm = TRUE)
    if(min_temp < min_abs) {min_abs <- min_temp}
  }
  return(min_abs)
}

plot_growth_curves <- function(dt, ids, col=rep("blue", length(ids)), title=NULL) {
  c = 1
  plt <- ggplot() + theme_classic() + ggtitle(title) + xlab("Time (h)") + ylab(bquote('Absorbance ('*OD[600]*')'))
  t_max <- nrow(dt)/max(dt$id)
  dt <- copy(dt); dt <- dt[, newtime := time*0.427718]
  
  for (i in ids) {
    plt <- plt + 
      geom_point(mapping = aes_string(x = dt[id == i]$newtime, y = dt[id == i]$abs), color=col[c]) + 
      scale_x_continuous(breaks=seq(0,72,12))
    #geom_line(mapping = aes_string(x = 1:t_max, y = dt[id == i]$fit)) 
    c <- c + 1
  }
  print(plt)
  return(plt)
}

plot_2microplates <- function(mp1, mp2, col=c("green","red")) {
  plt <- ggplot(data = mp1$data) +
    geom_point(mapping=aes(x=time, y=abs), color=col[1]) +
    geom_line(mapping=aes(x=time, y=fit)) +
    facet_wrap(~well, nrow=8)
  plt <- plt +
    geom_point(data=mp2$data, mapping=aes(x=time, y=abs), color=col[2]) +
    geom_line(data=mp2$data, mapping=aes(x=time, y=fit)) +
    facet_wrap(~well, nrow=8)
  print(plt)
}


process_microplate <- function(input_file) {
  # format raw plate reader (BMG SPECTROstar Omega) kinetic curve data using "gsed" command from the Linux terminal
  output_file <- paste0(strsplit(input_file, ".", fixed=TRUE)[[1]][1], "-cleaned.CSV")
  system( paste0("gsed -n '/[ABCDEFGH],/p' ", input_file, " >", output_file ) )

  # reads formatted .csv data
  df <- fread(paste0(script.dir, "/", output_file), header=FALSE)

  N_CYCLES <- dim(df)[1] / 8 #extract number of cycles

  row.names_ <- as.character(df[[1]][seq(1,8)])
  column.names_ <- as.character(sprintf('%02d',seq(1,12)))
  colnames(df) <- c("row", column.names_)

  df.melt <- melt(df, id.vars <- "row", value.name = "abs")
  colnames(df.melt)[[2]] <- "col"

  # assigns a unique ID to each well and reorder columns
  df.melt$rowIndex <- as.numeric(lapply(as.character(df.melt$row), utf8ToInt)) - 65
  df.melt$id <- as.numeric(df.melt$col) + (12 * df.melt$rowIndex)
  df.melt$well <- paste0(df.melt$row, df.melt$col)
  df.melt[, c("rowIndex", "row", "col") := NULL] # remove helper columns
  df.melt[, "time" := 1:.N , by = id] # create a time column (reverse order: .N:1)
  setcolorder(df.melt, c("id", "well", "time", "abs"))

  # creates a list (nls.list) of nonlinear least squares models with parametrized logistic curve of the form: " a/( 1 + e^(-b(t-t0)) ) + c "
  # because of warnOnly=TRUE, fitter will not throw errors even if model fitting failed. use with caution!
  nls.list <- vector("list", max(df.melt$id))
  if ( !(N_CYCLES<=2) ) {
    for (i in 1:max(df.melt$id)) {
      gcurve <- df.melt[id==i]
      nls.list[[i]] <- nls(abs ~ a/(1 + exp(-b*(time-t0)))+c, start=list(a=1, b=0.05, c=1, t0=N_CYCLES/2), data=gcurve, algorithm="port", control=nls.control(maxiter=1000, warnOnly=TRUE), lower=c(-3,-1,-10,-200), upper=c(3,1,10,700) )
    }

    df.melt[, "fit" := fit_nls_curve(nls.list[[id]], data.frame(time=1:.N), anyNA(abs)), by=id]

    return( new_microplate(df.melt, nls.list) )
  }
  else {
    print(paste0("Plate ", input_file, " is endpoint"))
    df.melt[, "fit" := abs, by=id]
    return( new_microplate(df.melt, endpoint = TRUE) )
  }
}




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


