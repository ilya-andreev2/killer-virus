library(DescTools)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(ggpubr)
library(gggenes)
library(gdata)

# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd() # is used in code - need to Source script, not Run
source("../utils.R")

# FUNCTIONS
prob_median <- function(x) {
  if (length(x) %% 2 == 1) {
    return(median(x))
  }
  else {
    return(x[sample(c(length(x) / 2, length(x) / 2 + 1), 1)])
  }
}

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

microplates <- list()
for (i in 1:length(file_names)) {
  microplates[[i]] <- process_microplate(file_names[[i]])
}

s <- data.table(
  plate = rep(1:20, each = 96),
  plate_type = rep(c("NT", "T"), each = 96 * 10),
  id = rep(1:96, 10),
  strain = rep(paste_expand("A", sprintf("%02d", 1:10), "_", sprintf("%02d", 1:96)), 2)
)
setkey(s, plate_type, plate, id, strain)

largeDT <- s[rep(seq_len(nrow(s)), each = 37), ]
largeDT[, "time" := rep(1:37, 20 * 96) * 1.6927927928]
names(largeDT)[3] <- "wid"
largeDT[, "abs" := microplates[[plate]]$data[id == wid, abs], by = .(plate, wid)] # without specifying 'by', variable (plate) will return a vector, not a single number -> assignment will fail
largeDT[, "abs" := abs - quantile(largeDT[, .SD[1:5], by = .(plate, wid)]$abs, 0.1)] # experiment-wide blank correction: quantile(, 0.5) would simply be the mean
setcolorder(largeDT, c("plate", "wid", "plate_type", "strain", "time", "abs"))

s_stat <- copy(s)
s_stat <- s_stat[plate <= 10]

pheno_fAUC <- largeDT[, AUC(time, abs), by = .(plate, strain)]
pheno_fAUC[strain %in% pheno_fAUC[plate <= 10 & V1 < 10]$strain, V1 := NA] # exclude poorly-growing strains
pheno_fAUC <- pheno_fAUC[, lapply(.SD, function(x) {
  return(x[2] / x[1])
}), by = .(strain), .SDcols = "V1"]
s_stat[, "pheno_fAUC" := pheno_fAUC[, 2]]

# reformatting the data for QTL mapping
pheno.K28 <- as.data.frame(s_stat$pheno_fAUC)
rownames(pheno.K28) <- s_stat$strain
pheno.K28 <- t(pheno.K28)

get.LOD.by.COR <- function(n.pheno, pheno, gdata) {
  # Lynch and Walsh p. 454
  return((-n.pheno * log(1 - cor(pheno, gdata, use = "pairwise.complete.obs")^2)) / (2 * log(10)))
}

# load segregant genotype data: -1  = BY, +1  = RM
load("1000BYxRM_with_names.RData")

# match genotype with phenotype: pheno_input_matrix = pheno.data -- geno_matrix = BYxRM_orig
match_pheno_and_geno <- function(pheno_input_matrix, geno_matrix) {
  BYxRM_strain_name <- do.call("rbind", strsplit(rownames(geno_matrix), ":"))[, 1]
  sname <- (do.call("rbind", strsplit(colnames(pheno_input_matrix), "-"))[, 1])
  pheno_input_matrix <- pheno_input_matrix[, -which(is.na(match(sname, BYxRM_strain_name)))]
  sname <- (do.call("rbind", strsplit(names(pheno_input_matrix), "-"))[, 1])
  gdata <- geno_matrix[match(sname, BYxRM_strain_name), ]
  return(list(pheno = pheno_input_matrix, gdata = gdata))
}

pg <- match_pheno_and_geno(pheno.K28, BYxRM_orig)
p <- t(pg$pheno)

# removes non-informative markers
g <- pg$gdata
g_dd <- g[, -which(duplicated(g, MARGIN = 2))] # dd for 'deduplicated' - SNPs with duplicate genotypes removed
pp <- apply(p, 2, as.numeric)
N <- length(pp)

#### QTL mapping ####
LODs <- get.LOD.by.COR(N, pp, g)

# histogram / distribution of fAUC phenotype
hist(pheno.K28)

# Confidence interval by bootstrap with n=1000 samplings
# !uncomment the code below if you want to run it yourself
#
# nboot <- 1000 # number of bootstrap samples
# peak_vector <- vector("character")
# for (i in 1:nboot) { # bootstrap QTL mapping, N times, sample strains with replacement
#   g_index <- base::sample(1:N, replace = TRUE)
#   g_s = g[g_index,]; pp_s <- pp[g_index] # g with duplications of genotypes or g_dd (without)
#   LODs_s <- get.LOD.by.COR(N, pp_s, g_s)
#   peak_index <- prob_median(which(LODs_s == max(LODs_s)))
#   peak_vector <- c(peak_vector, names(LODs_s[, peak_index]))
#   print(i)
# }
# tmp <- sapply(peak_vector, parse_coord)
# table(tmp)

# quantile(tmp, c(0.025, 0.975))
d_max025 <- 180697 # d_max025 = quantile(tmp, .025);
d_max975 <- 184686 # d_max975 = quantile(tmp, .975);

CI_bootstrap <- c(d_max025, d_max975)
peak_coords <- tstrsplit(names(LODs[, which.max(LODs)]), "_")[3] # chrI peak coordinate
print(sprintf("chrI peak coordinate(s): %.0f", peak_coords))
print(sprintf("Confidence interval (bootstrap): [%.0f, %.0f]", CI_bootstrap[1], CI_bootstrap[2]))

# extract genotypes for QTLs with max LOD score
max_ind <- which.max(LODs)
match_ind <- match(sapply(rownames(g), function(x) {
  return(tstrsplit(x, split = ":")[[1]])
}), colnames(pheno.K28))

peak_geno_pheno <- data.table(
  geno = as.factor(ifelse(g[, max_ind] == -1, "BY", "RM")),
  pheno = pheno.K28[, match_ind]
)

# Compare group means with Welch 2-sample t-test
t.test(peak_geno_pheno[geno == "BY", pheno], peak_geno_pheno[geno == "RM", pheno])

# Using ggpubr
(p <- ggboxplot(peak_geno_pheno, x = "geno", y = "pheno", add = "jitter", seed = 101010) + # what exactly do the low/high bars represent?
  stat_compare_means(comparisons = list(c("BY", "RM")), label.y = 1.4) +
  labs(x = "Genotype at chrI locus", y = "fAUC Phenotype"))
ggsave("boxplot_jitter.pdf", device = "pdf", width = 5, height = 6.3)

# peak_geno_pheno[, quantile(pheno, 0.75, na.rm = TRUE), by=.(geno)] # to determine what's considered "dead"


#### Manhattan plot with ggplot2 ####
x <- unlist(tstrsplit(colnames(g), split = "_")[2])
LODs.dt <- as.data.table(t(LODs), keep.rownames = TRUE)
LODs.dt[, "g_pos" := parse_coord(rn)][, "chr" := ..x]
colnames(LODs.dt)[1:3] <- c("variant", "LOD", "g_pos")
breaks <- parse_coord(c(colnames(g)[-which(duplicated(x))], last(colnames(g_dd))))
len_b <- length(breaks)
mid_breaks <- vector("numeric", len_b - 1)
for (i in 1:(len_b - 1)) {
  mid_breaks[i] <- (breaks[i] + breaks[i + 1]) / 2
}

sig_thr <- 3.75 # 1000 permutations significance threshold
#### !1000 permutations threshold code ####
# permu=matrix(NA, nrow = 1000, ncol=2)
# for (j in 1:1000){
#   print(j)
#   this_p = pp[sample(1:length(pp), replace=TRUE)]
#   this_n = length(this_p)
#   middle = get.LOD.by.COR(this_n, this_p, g)
#   permu[j,]=apply(middle, 1, max, na.rm=TRUE)
# }
# #
# maxquan=apply(permu,2,quantile,0.95)

# Whole-genome Manhattan Plot
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
ggsave("Manhattan_Plot.pdf", device = "pdf", width = 6, height = 3.5) # export
