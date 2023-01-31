library(DescTools)
library(ggpubr)
library(gdata)

# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")

# FUNCTIONS

# prob_median: returns the median element if the input vector is odd length, or randomly one randomly selected out of two median elements if the length is even
prob_median <- function(x) {
  if (length(x) %% 2 == 1) {
    return(median(x))
  }
  else {
    return(x[sample(c(length(x) / 2, length(x) / 2 + 1), 1)])
  }
}

# parse_coord: parses genomic coordinates from name
parse_coord <- function(x) {
  return(as.double(unlist(tstrsplit(x, split = "_")[1])))
}

# MAIN SCRIPT
file_names <- list(
  "A1_NT_02-14-2020.CSV", "A2_NT_02-14-2020.CSV",
  "A3_NT_02-14-2020.CSV", "A4_NT_02-14-2020.CSV",
  "A5_NT_02-14-2020.CSV", "A6_NT_02-14-2020.CSV",
  "A7_NT_02-14-2020.CSV", "A8_NT_02-14-2020.CSV",
  "A9_NT_02-14-2020.CSV", "A10_NT_02-14-2020.CSV",
  "A1_T_02-14-2020.CSV", "A2_T_02-14-2020.CSV",
  "A3_T_02-14-2020.CSV", "A4_T_02-14-2020.CSV",
  "A5_T_02-14-2020.CSV", "A6_T_02-14-2020.CSV",
  "A7_T_02-14-2020.CSV", "A8_T_02-14-2020.CSV",
  "A9_T_02-14-2020.CSV", "A10_T_02-14-2020.CSV"
)

microplates <- list() #initialize empty plate list
for (i in 1:length(file_names)) { #process 96-well plate data from BMG SpectroSTAR Omega using a custom function defined in 'utils.R'
  microplates[[i]] <- process_microplate(file_names[[i]])
}


s <- data.table( #define a data table, similar to a data frame but more powerful
  plate = rep(1:20, each = 96), # there are 20 plates - see 'file_names" above
  plate_type = rep(c("NT", "T"), each = 96 * 10), #half are with (T)oxin, the other half are (N)o-(T)oxin
  id = rep(1:96, 10), #each plate has 96 wells
  strain = rep(paste_expand("A", sprintf("%02d", 1:10), "_", sprintf("%02d", 1:96)), 2) 
)
setkey(s, plate_type, plate, id, strain) #order by plate_type, then by plate, then by id, then by strain

largeDT <- s[rep(seq_len(nrow(s)), each = 37), ] #replicate each entry 37 times because there are 37 time points
largeDT[, "time" := rep(1:37, 20 * 96) * 1.6927927928] #convert index to real time: take the very last endpoint (62h38m, or 62.63h) in the last read plate, then the conversion factor is 62.63 h / 37 = 1.6927927928
names(largeDT)[3] <- "wid" #well ID
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)] # without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, "abs" := abs - quantile(largeDT[, .SD[1:5], by = .(plate, wid)]$abs, 0.1)] # experiment-wide blank correction: quantile(, 0.5) would simply be the mean
setcolorder(largeDT, c("plate", "wid", "plate_type", "strain", "time", "abs"))
largeDT #now we have normalized absorbance over absolute time - can do analysis!

# Compute fAUC (=AUC_+K28 / AUC_-K28) phenotype
pheno_fAUC <- largeDT[, AUC(time, abs), by = .(plate, strain)] #compute AUC using trapezoid method
pheno_fAUC[strain %in% pheno_fAUC[plate <= 10 & V1 < 10]$strain, V1 := NA] # identify poorly-growing strains
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x) {
  return(x[2] / x[1])
}), by = .(strain), .SDcols = "V1"] #groups by strain (two rows in each group: NT and T respectively), then divides AUC_T by AUC_NT

s_stat <- copy(s) # another data.table skeleton
s_stat <- s_stat[plate <= 10]
s_stat[, "pheno_fAUC" := pheno_fAUC[,2]] # assigns a new column "pheno_fAUC" and copies values
s_stat <- s_stat[!is.na(pheno_fAUC)] # exclude 34 poorly growing strains from analysis: 960 -> 926 (check: s_stat[, .N])

# reformatting the data for QTL mapping algorithm
pheno.K28 <- as.data.frame(s_stat$pheno_fAUC)
rownames(pheno.K28) <- s_stat$strain
pheno.K28 <- t(pheno.K28)

# A portion of the following code was reused with permission from Gou et al 'The Genetic Basis of Mutation Rate Variation in Yeast' (Genetics, 2019)
# Original data is also available at https://github.com/gouliangke/Mutation-rate/tree/master/genotype)
get.LOD.by.COR <- function(n.pheno, pheno, gdata) {
  # Lynch and Walsh p. 454
  return((-n.pheno * log(1 - cor(pheno, gdata, use = "pairwise.complete.obs")^2)) / (2 * log(10)))
}

# load segregant genotype data: -1  = BY, +1  = RM
load("1000BYxRM_with_names.RData")

# match genotype with phenotype: pheno_input_matrix = pheno.data; geno_matrix = BYxRM_orig
match_pheno_and_geno <- function(pheno_input_matrix, geno_matrix) {
  BYxRM_strain_name <- do.call("rbind", strsplit(rownames(geno_matrix), ":"))[, 1]
  sname <- (do.call("rbind", strsplit(colnames(pheno_input_matrix), "-"))[, 1])
  pheno_input_matrix <- pheno_input_matrix[, -which(is.na(match(sname, BYxRM_strain_name)))]
  sname <- (do.call("rbind", strsplit(names(pheno_input_matrix), "-"))[, 1])
  gdata <- geno_matrix[match(sname, BYxRM_strain_name), ]
  return(list(pheno = pheno_input_matrix, gdata = gdata))
}

pg <- match_pheno_and_geno(pheno.K28, BYxRM_orig) # some strains have no match - 926 -> 912 strains at this step. (check: length(p))
p <- t(pg$pheno)

# remove non-informative markers
g <- pg$gdata
g_dd <- g[, -which(duplicated(g, MARGIN = 2))] # dd for 'deduplicated' - SNPs with duplicate genotypes are removed
pp <- apply(p, 2, as.numeric)
N <- length(pp)

#### QTL mapping ####
LODs <- get.LOD.by.COR(N, pp, g)
LODs_dd <- get.LOD.by.COR(N, pp, g_dd) # the results are the same whether you use g or g_dd

LODs_dt <- as.data.table(t(LODs), keep.rownames = TRUE)
LODs_dt[, i := 1:length(LODs)]
LODs_dt[, gpos := parse_coord(rn)]
LODs_dt[, chr := unlist(tstrsplit(rn, split = "_")[2])]
LODs_dt[, chrpos := as.double(unlist(tstrsplit(rn, split = "_")[3]))]
LODs_dt[chr == "chrI"][prob_median(which(V1 == max(V1)))]

# histogram / distribution of fAUC phenotype
hist(pheno.K28)

# Confidence interval by bootstrap with n=1000 samplings
nboot <- 1000 # number of bootstrap samples
peak_vector_I <- vector("integer")
peak_vector_XII <- vector("integer")
peak_vector_XIV <- vector("integer")
for (i in 1:nboot) { # bootstrap QTL mapping, N times, sample strains with replacement
  g_index <- base::sample(1:N, replace = TRUE)
  g_s = g[g_index,]; pp_s <- pp[g_index] # g with duplications of genotypes or g_dd (without)
  LODs_s <- get.LOD.by.COR(N, pp_s, g_s)
  LODs_s_dt <- as.data.table(t(LODs_s), keep.rownames = TRUE)
  LODs_s_dt[, i := 1:length(LODs_s)]
  LODs_s_dt[, gpos := parse_coord(rn)]
  LODs_s_dt[, chr := unlist(tstrsplit(rn, split = "_")[2])]
  LODs_s_dt[, chrpos := as.double(unlist(tstrsplit(rn, split = "_")[3]))]
  peak_vector_I <- c(peak_vector_I, LODs_s_dt[chr == "chrI"][prob_median(which(V1 == max(V1))), chrpos])
  peak_vector_XII <- c(peak_vector_XII, LODs_s_dt[chr == "chrXII"][prob_median(which(V1 == max(V1))), chrpos])
  peak_vector_XIV <- c(peak_vector_XIV, LODs_s_dt[chr == "chrXIV"][prob_median(which(V1 == max(V1))), chrpos])
  print(i)
}
tmp <- sapply(peak_vector_I, parse_coord)
table(tmp)

# Compute 95% confidence intervals (for Supplementary Table S2)
dI_max025 = quantile(peak_vector_I, .025, type = 1)
dI_max975 = quantile(peak_vector_I, .975, type = 1)
dXII_max025 = quantile(peak_vector_XII, .025, type = 1)
dXII_max975 = quantile(peak_vector_XII, .975, type = 1)
dXIV_max025 = quantile(peak_vector_XIV, .025, type = 1)
dXIV_max975 = quantile(peak_vector_XIV, .975, type = 1)

LODs_dt[, .SD[V1 == max(V1)], by="chr"] # peak markers by chromosome

CI_bootstrap <- c(dI_max025, dI_max975) #chrI peak 95% CI: 180,646 to 184,686
peak_coords <- tstrsplit(names(LODs[, which.max(LODs)]), "_")[3] # chrI peak coordinate
print(sprintf("chrI peak coordinate(s): %.0f", peak_coords))
print(sprintf("Confidence interval (bootstrap): [%.0f, %.0f]", CI_bootstrap[1], CI_bootstrap[2]))
print(sprintf("Confidence interval length: %.0f bp", CI_bootstrap[2] - CI_bootstrap[1]))

# Extract genotypes for QTLs with max LOD score
max_ind <- which.max(LODs) #tie-break: which.max takes the first max occurrence
match_ind <- match(sapply(rownames(g), function(x) {
  return(tstrsplit(x, split = ":")[[1]])
}), colnames(pheno.K28))

peak_geno_pheno <- data.table(
  geno = as.factor(ifelse(g[, max_ind] == -1, "BY", "RM")),
  pheno = pheno.K28[, match_ind]
)
peak_geno_pheno <- peak_geno_pheno[!is.na(pheno)] # remove NA's

# Compare group means using Welch 2-sample t-test
t.test(peak_geno_pheno[geno == "BY", pheno], peak_geno_pheno[geno == "RM", pheno])

# % segregants alive, grouped by genotype, fAUC = 0.5 cutoff
peak_geno_pheno[, .SD[pheno > 0.5, .N]/.N, by = geno]

# Pre-processing for Manhattan plot
x <- unlist(tstrsplit(colnames(g), split = "_")[2])
LODs.dt <- as.data.table(t(LODs), keep.rownames = TRUE)
LODs.dt[, "g_pos" := parse_coord(rn)][, "chr" := ..x]
colnames(LODs.dt)[1:3] <- c("variant", "LOD", "g_pos")
breaks <- parse_coord(c(colnames(g)[-which(duplicated(x))], last(colnames(g_dd))))
len_b <- length(breaks)
mid_breaks <- vector("numeric", len_b - 1)
for (i in 1:(len_b - 1)) {
  mid_breaks[i] <- (breaks[i] + breaks[i+1]) / 2
}

# Calculating 5% FWER via 1000 permutations test
# permu=matrix(NA, nrow = 1000, ncol=2)
# for (j in 1:1000){
#   print(j)
#   this_p = pp[sample(1:length(pp))]
#   this_n = length(this_p)
#   middle = get.LOD.by.COR(this_n, this_p, g)
#   permu[j,]=apply(middle, 1, max, na.rm=TRUE)
# }
# maxquan = sort(permu[,2])[950]
sig_thr <- 3.552774 # pre-defined as a "magic number" from one of the code runs - uncomment code above to compute your own
#sig_thr <- maxquan # if the code above was run

# Figure 2B: Whole-genome Manhattan Plot (LOD plot)
(manhplot <- ggplot() +
  geom_path(LODs.dt, mapping = aes(x = g_pos, y = LOD), alpha = 0.85, color = "blue", size = 0.3) +
  geom_hline(yintercept = sig_thr, color = "grey", alpha = 0.7) +
  geom_text(aes(breaks[len_b - 1], sig_thr, label = "5% FWER", vjust = -0.75), size = 2.5, color = "grey") +
  scale_x_continuous(
    expand = c(0, 0),
    label = c(rep("", len_b), tstrsplit(unique(LODs.dt$chr), "chr")[[2]]),
    breaks = c(breaks, mid_breaks)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, max(LODs) + 5)) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = "Genome Position (chr)", y = "LOD Score") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8),
    axis.ticks.x = element_line(color = c(rep("black", len_b), rep(NA, len_b)))
  )
)
ggsave("T10_LODplot.pdf", device = "pdf", width = 6, height = 4.6)


# Figure 2C
(p <- ggboxplot(peak_geno_pheno, x = "geno", y = "pheno", add = "jitter", seed = 101010) + 
   stat_compare_means(comparisons = list(c("BY", "RM")), label.y = 1.4) +
   labs(x = "Genotype at chrI locus", y = "fAUC Phenotype"))
ggsave("T10_pheno_boxplot.pdf", device = "pdf", width = 5, height = 6.3)
