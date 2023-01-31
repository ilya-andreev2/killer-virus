# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(gdata)
library(stringr)
library(data.table)

# OBJECTS AND METHODS
new_microplate <- function(data = data.table(),
                           nls_list = list(),
                           endpoint = FALSE) {
  # <<<<<<<<<TO IMPLEMENT: name, date, !cycle_time!, protocol, other metadata from .csv>>>>>>>>>>
  stopifnot(is.data.table(data))
  plate <- list(data = data, nls_list = nls_list, endpoint = endpoint)
  attr(plate, "class") <- "microplate"
  plate
}

print.microplate <- function(microplate, col = "magenta", ...) {
  plt <- ggplot(data = microplate$data) +
    geom_point(mapping = aes(x = time, y = abs), color = col) +
    facet_wrap(~ well, nrow = 8) #+
  print(plt)
  print(microplate$data)
}

# FUNCTIONS
blank_correct <- function(s) {
  min_abs <- 999999
  for (i in seq_len(nrow(s))) {
    print(min(microplates[[s[i, plate]]]$data[id == s[i, id]]$abs))
    min_temp <- min(microplates[[s[i, plate]]]$data[id == s[i, id]]$abs, na.rm = TRUE)
    if (min_temp < min_abs) {
      min_abs <- min_temp
    }
  }
  return(min_abs)
}

plot_growth_curves <- function(dt, ids, time_factor = 0.427718, col = rep("blue", length(ids)), title = NULL) {
  c <- 1
  plt <- ggplot() +
    theme_classic() +
    ggtitle(title) +
    xlab("Time (h)") +
    ylab(bquote("Absorbance (" * OD[600] * ")"))
  dt <- copy(dt)
  dt <- dt[, newtime := time * time_factor]
  for (i in ids) {
    plt <- plt +
      geom_point(mapping = aes_string(x = dt[id == i]$newtime, y = dt[id == i]$abs), color = col[c]) +
      scale_x_continuous(breaks = seq(0, 72, 12))
    c <- c + 1
  }
  print(plt)
  return(plt)
}

plot_2microplates <- function(mp1, mp2, col = c("green", "red")) {
  plt <- ggplot(data = mp1$data) +
    geom_point(mapping = aes(x = time, y = abs), color = col[1]) +
    facet_wrap(~well, nrow = 8)
  plt <- plt +
    geom_point(data = mp2$data, mapping = aes(x = time, y = abs), color = col[2]) +
    facet_wrap(~well, nrow = 8)
  print(plt)
}

process_microplate <- function(input_file) {
  # format raw plate reader (BMG SPECTROstar Omega) kinetic curve data using "gsed" command from the Linux terminal
  output_file <- paste0(strsplit(input_file, ".", fixed = TRUE)[[1]][1], "-cleaned.CSV")
  system(paste0("gsed -n '/[ABCDEFGH],/p' ", input_file, " >", output_file)) #requires sed installed. For Win users, use sed from GnuWin32 and replace single quotes with double quotes in regex

  # reads formatted .csv data
  df <- fread(paste0(script_dir, "/", output_file), header = FALSE)
  column.names_ <- as.character(sprintf("%02d", seq(1, 12)))
  colnames(df) <- c("row", column.names_)

  df.melt <- melt(df, id.vars = "row", value.name = "abs")
  colnames(df.melt)[[2]] <- "col"

  # assigns a unique ID to each well and reorder columns
  df.melt$rowIndex <- as.numeric(lapply(as.character(df.melt$row), utf8ToInt)) - 65
  df.melt$id <- as.numeric(df.melt$col) + (12 * df.melt$rowIndex)
  df.melt$well <- paste0(df.melt$row, df.melt$col)
  df.melt[, c("rowIndex", "row", "col") := NULL] # remove helper columns
  df.melt[, "time" := 1:.N, by = id] # create a time column (reverse order: .N:1)
  setcolorder(df.melt, c("id", "well", "time", "abs"))

  return(new_microplate(df.melt))
}

paste_expand <- function(...) {
  args_ <- list(...)
  len <- lengths(args_)
  N <- prod(len)
  result <- vector("character", N)
  for (i in seq_len(length(args_))) {
    arg <- as.character(args_[[i]])
    N <- N / len[i]
    result <- paste0(result, rep(arg, each = N))
  }
  return(result)
}
