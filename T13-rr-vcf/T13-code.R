library(data.table)
library(VariantAnnotation)
library(GenomicFeatures)
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd() # is used in code - need to Source script, not Run
options(showHeadLines=1200); options(max.print = 999999); options(width = 10000)
#options(max.print = 100)

# FUNCTIONS
process_variants <- function(vcf, cdb, strains, genes) {
  vcf.geno <- geno(vcf)$GT
  var <- cdb[ which(cdb$GENEID %in% genes) ]
  varnames <- names(var)
  for (s in strains) {
    elementMetadata(var)[[ as.character(s) ]] <- vcf.geno[varnames,s] #adds the genotype columns from vcf to var
  }
  var.df <- data.frame(var[,c("GENEID","PROTEINLOC","CONSEQUENCE","REFCODON","VARCODON","REFAA","VARAA",strains)])
  var.df$names <- names(var); var.df <-var.df[,c("names",colnames(var.df)[1:(length(colnames(var.df))-1)])]
  return(var.df)
}


# MAIN SCRIPT
vcf <- readVcf("RR_parents.vcf") # this file can be obtained from Joshua Bloom (jbloom@mednet.ucla.edu); citation: doi.org/10.7554/eLife.49212
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
coding <- predictCoding(vcf, txdb, seqSource=Scerevisiae)


# extract variants within YAR027/28W for BY and RM
gene_list <- c("YAR027W","YAR028W")
strains <- c("BYa","RMx")
variants_BYxRM.df <- process_variants(vcf, coding, strains, gene_list)

# TO-DO: Fix YAR023C ranges in the vcf file
# select(txdb, keys = "YAR023C", columns=columns(txdb), keytype="GENEID")
# txdb_dump <- as.list(txdb)
# txdb_dump[[1]][which(txdb_dump[[1]]$tx_name == "YAR023C"),]
# new_txdb <- do.call(makeTxDb, txdb_dump)
# Not gonna work - more upstream start codon at pos 180,094 (Ref), which is where RM starts its transcript, is out of frame with the annotated ORF in Ref genome - would have to reconcile it with the .vcf file and find where the frameshift had occurred, etc.


# extending the variant table to all DUP240 genes
DUP240_list <- c("YAR023C","YAR027W","YAR028W","YAR029W","YAR031W","YAR033W","YCR007C","YGL053W","YGL051W","YHL044W")
strains <- c("M22","BYa","RMx","YPS163a","YJM145x","CLIB413a","YJM978x","YJM454a","YPS1009x","I14a","Y10x","PW5a","273614xa","YJM981x","CBS2888a","CLIB219x")
variants_DUP240.df <- process_variants(vcf, coding, strains, DUP240_list)
variants_DUP240.df[ which(variants_DUP240.df$CONSEQUENCE != "synonymous"), ]

# Print out your .df in console, copy-paste output into a .txt, then change to the same font in your text editor to have all columns line up nicely
#write.table( x = variants.df, file = "variants.df.csv", sep=";", col.names=TRUE, row.names=FALSE, quote=FALSE )

# peak 2 from BYxRM QTL mapping
peak2_gene_list <- c("YNL087W","YNL086W","YNL085W","YNL084C","YNL083W","YNL082W","YNL081C","YNL080C","YNL079C","YNL078W","YNL077W")
strains <- c("BYa","RMx")
variants_peak2.df <- process_variants(vcf, coding, strains, peak2_gene_list)
variants_peak2.df[ which(variants_peak2.df$CONSEQUENCE != "synonymous"), ]

# YJM454xYPS1009
MATlocus_gene_list <- c("YCR028C","YCR030C","YCR031C","YCR032W","YCR033W","YCR034W","YCR035C","YCR036W","YCR037C","YCR038C")
strains <- c("BYa","YJM454a","YPS1009x")
variants_MAT.df <- process_variants(vcf, coding, strains, MATlocus_gene_list)
variants_MAT.df[ which((variants_MAT.df$CONSEQUENCE != "synonymous")), ]

#for Figure 3
strains <- c("M22","BYa","RMx","YPS163a","YJM145x","CLIB413a","YJM978x","YJM454a","YPS1009x","I14a","Y10x","PW5a","273614xa","YJM981x","CBS2888a","CLIB219x")
vKTD1.df <- process_variants(vcf, coding, strains, "YAR028W")
vKTD1.df$PROTEINLOC <- as.integer(unlist(vKTD1.df[,8]))

vKTD1.melt <- as.data.table(reshape2::melt(vKTD1.df, id.vars = c("PROTEINLOC", "CONSEQUENCE"), measure.vars = colnames(vKTD1.df)[14:(14+15)]))
vKTD1.melt <- vKTD1.melt[!(variable %in% c("M22", "YJM145x", "YJM454a", "PW5a", "YJM978x", "CLIB219x"))]
vKTD1.melt[value == "2", "value" := "1"] # set alt-2 alleles to 1
vKTD1.melt[, "value" := as.integer(value)]
vKTD1.melt$PROTEINLOC <- unlist(vKTD1.melt$PROTEINLOC)
vKTD1.melt <- vKTD1.melt[value > 0 | variable == "BYa"]

# PROBLEM: Missing variant at AA 12 in YPS1009/CBS2888
vKTD1.melt <- rbind(vKTD1.melt, 
      data.table(PROTEINLOC=c(12,12), CONSEQUENCE=rep("nonsynonymous", 2), variable=c("YPS1009x", "CBS2888a"), value=c(1,1))
)

unique(vKTD1.melt$variable)

g <- ggplot(data = vKTD1.melt, mapping = aes(x = PROTEINLOC, y = value, fill = CONSEQUENCE)) + 
   scale_x_continuous(breaks = c(seq(0,234,50), 234)) + 
   #geom_point() +
   geom_rect(aes(xmin=PROTEINLOC, xmax=PROTEINLOC+0.25, ymin=0.95, ymax=1.05)) +
   facet_grid(rows = factor(vKTD1.melt$variable, levels=c("X273614xa", "Y10x", "CLIB413a", "YPS163a", "BYa", "I14a", "YJM981x", "RMx", "CBS2888a", "YPS1009x", "M22", "YJM145x", "YJM978x", "YJM454a", "PW5a", "CLIB219x"))) + 
   coord_cartesian(
      xlim = c(0,234),
      ylim = c(0.9, 1.1)
   ) +
   theme_minimal() +
   theme(axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_line(),
         axis.line.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(fill = NA),
         strip.text.y = element_text(angle = 0)
   )
g

# check theme properties of ggplot
g$theme

ggsave("Fig3A_variants.pdf", width = 6, height = 9)

# Analysis of Susheela's genes mutations in RR strains
res_list <- read.csv("carroll_resistant_muts.csv", header=FALSE, fileEncoding="UTF-8-BOM")
res_list <- res_list$V1
strains <- c("M22","BYa","RMx","YPS163a","YJM145x","CLIB413a","YJM978x","YJM454a","YPS1009x","I14a","Y10x","PW5a","273614xa","YJM981x","CBS2888a","CLIB219x")
variants_res.df <- process_variants(vcf, coding, strains, res_list)
variants_res.df <- variants_res.df[ which(variants_res.df$CONSEQUENCE != "synonymous"), ]
variants_res.df <- variants_res.df[ which(variants_res.df$CONSEQUENCE != "nonsynonymous"), ]
#fwrite(variants_res.df, file="all_res_muts.txt", sep="\t")
View(variants_res.df)

sens_list <- read.csv("carroll_sensitive_muts.csv", header=FALSE, fileEncoding="UTF-8-BOM")
sens_list <- sens_list$V1
strains <- c("M22","BYa","RMx","YPS163a","YJM145x","CLIB413a","YJM978x","YJM454a","YPS1009x","I14a","Y10x","PW5a","273614xa","YJM981x","CBS2888a","CLIB219x")
variants_sens.df <- as.data.table(process_variants(vcf, coding, strains, sens_list))
variants_sens.df <- variants_sens.df[ which(variants_sens.df$CONSEQUENCE != "synonymous"), ]
variants_sens.df <- variants_sens.df[ which(variants_sens.df$CONSEQUENCE != "nonsynonymous"), ]
#fwrite(variants_sens.df, file="all_sens_muts.txt", sep="\t")
View(variants_sens.df)

# Manually curated for "fake" frameshifts, verified by Pacbio sequence analysis (Meru's algo)
# some Meru-confirmed are b/c they're only a few aa shorter than reference full-length
# YBL041W - only 5 aa shorter
# YER151C - fake 1st codon frameshifts
# YFL023W - co-occurrence + 9,6 lengths
# YGL008C - co-occurrence + 18,6 lengths
# YIL110W - fake last codon variant
# YJR093C - Meru confirmed OK
# YKL114C - Meru confirmed OK
# YLR175W - Meru confirmed OK
# YOR262W - only 5 amino acids shorter
# YOR310C - Meru confirmed OK
# YPR133C - co-occurrence + 9,9 lengths
variants_sens.df2 <- variants_sens.df[!(GENEID %in% c("YBL041W", "YER151C", "YFL023W", "YGL008C", "YIL110W", "YJR093C", "YKL114C", "YLR175W", "YOR262W", "YOR310C", "YPR133C"))]
fwrite(variants_sens.df2, file="all_sens_muts_edited.txt", sep="\t")
View(variants_sens.df2)

variants_sens.df3 <- variants_sens.df[GENEID %in% c("YBL041W", "YER151C", "YFL023W", "YGL008C", "YIL110W", "YJR093C", "YKL114C", "YLR175W", "YOR262W", "YOR310C", "YPR133C")]
