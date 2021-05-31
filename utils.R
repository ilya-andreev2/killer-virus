library(ggplot2)
library(gdata)
library(stringr)
library(data.table)

# OBJECTS
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
    #geom_line(mapping=aes(x=time, y=fit)) +
    facet_wrap(~well, nrow=8) #+
  #coord_cartesian(xlim=c(0,500))
  print(plt)
  print(microplate$data)
}

# FUNCTIONS #meow
# fit_nls_curve <- function(nls.model, time.points, na.flag) {
#   if(na.flag == TRUE) { return( as.double(rep(NA, length(time.points))) ) }
#   else { return( predict(nls.model, newdata=time.points) ) }
# }
# 
# refit_nls_curves <- function(p, i, a,b,c,t0) {
#   gcurve <- p$data[id==i]
#   p$nls_list[[i]] <- nls(abs ~ a/(1 + exp(-b*(time-t0)))+c, start=list(a=a, b=b, c=c, t0=t0), data=gcurve, algorithm="port", control=nls.control(maxiter=1000, warnOnly=TRUE), lower=c(-3,-1,-10,-200), upper=c(3,1,10,500) )
#   p$data[id==i, "fit" := fit_nls_curve(p$nls_list[[i]], data.frame(time=1:.N), anyNA(abs))]
#   return(p)
# }

blank_correct <- function(s) {
  min_abs <- 999999
  for (i in 1:nrow(s)) {
    print(min(microplates[[ s[i,plate] ]]$data[ id == s[i,id] ]$abs))
    min_temp <- min(microplates[[ s[i,plate] ]]$data[ id == s[i,id] ]$abs, na.rm = TRUE)
    if(min_temp < min_abs) {min_abs <- min_temp}
  }
  return(min_abs)
}

plot_growth_curves <- function(dt, ids, time_factor = 0.427718, col=rep("blue", length(ids)), title=NULL) {
  c = 1
  plt <- ggplot() + theme_classic() + ggtitle(title) + xlab("Time (h)") + ylab(bquote('Absorbance ('*OD[600]*')'))
  t_max <- nrow(dt)/max(dt$id)
  dt <- copy(dt); dt <- dt[, newtime := time*time_factor]
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
    #geom_line(mapping=aes(x=time, y=fit)) +
    facet_wrap(~well, nrow=8)
  plt <- plt +
    geom_point(data=mp2$data, mapping=aes(x=time, y=abs), color=col[2]) +
    #geom_line(data=mp2$data, mapping=aes(x=time, y=fit)) +
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
  # nls.list <- vector("list", max(df.melt$id))
  # if ( !(N_CYCLES<=2) ) {
  #   for (i in 1:max(df.melt$id)) {
  #     gcurve <- df.melt[id==i]
  #     nls.list[[i]] <- nls(abs ~ a/(1 + exp(-b*(time-t0)))+c, start=list(a=1, b=0.05, c=1, t0=N_CYCLES/2), data=gcurve, algorithm="port", control=nls.control(maxiter=1000, warnOnly=TRUE), lower=c(-3,-1,-10,-200), upper=c(3,1,10,700) )
  #   }
  #   df.melt[, "fit" := fit_nls_curve(nls.list[[id]], data.frame(time=1:.N), anyNA(abs)), by=id]
  #   return( new_microplate(df.melt, nls.list) )
  # }
  # else {
  #   print(paste0("Plate ", input_file, " is endpoint"))
  #   df.melt[, "fit" := abs, by=id]
  #   return( new_microplate(df.melt, endpoint = TRUE) )
  # }
  return( new_microplate(df.melt) )
}

paste_expand <- function(...) {
  args_ <- list(...)
  len <- lengths(args_)
  N <- prod(len)
  result <- vector("character", N)
  for (i in 1:length(args_)) {
    arg <- as.character(args_[[i]])
    N <- N / len[i]
    result <- paste0(result, rep(arg, each = N))
  }
  return(result)
}
