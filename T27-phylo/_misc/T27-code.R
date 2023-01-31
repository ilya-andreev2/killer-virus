library(readr)
library(data.table)
library(ape)
library(msa)
library(ggtree)
library(ggrepel)
library(treeio)

# sets working directory to the script location - must use RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#RR_DUP240 <- readDNAStringSet(paste0(getwd(),"/RR_DUP240/RR_DUP240.fsa")) #all alleles, no stop codon in KTD1
RR_DUP240 <- readDNAStringSet(paste0(getwd(),"/RR_DUP240/RR_DUP240_allAlleles.fsa"))
#RR_DUP240 <- subseq(RR_DUP240, end=width(RR_DUP240)-3) # remove stop codon
names(RR_DUP240) <- paste(names(RR_DUP240), " ")# extra space required for conversion to Phylip below
#RR_DUP240_prot <- translate(RR_DUP240)
writeXStringSet(RR_DUP240_prot, paste0(getwd(),"/RR_DUP240/RR_DUP240_p.fsa"))
nseq <- length(RR_DUP240) #number of DUP240s analyzed

# Multiple sequence alignment - may take a while
alignment_n <- DNAMultipleAlignment(msa(RR_DUP240, order="input"), use.names=TRUE) # convert to MultipleAlignment object
write.phylip(alignment_n, "RR_DUP240_allAlleles.aln.nuc")
alignment_p <- AAMultipleAlignment(msa(RR_DUP240_prot, verbose=TRUE, order="input"), use.names=TRUE)
write.phylip(alignment_p, "RR_DUP240.aln.prot")

#ape::write.dna(alignment_n, "RR_DUP240.aln.nuc", nbcol=1, colsep=" ", colw=1000000)


# Nucleotide-based tree using FastTree + visualization in ggtree (see treedata book)
#system("source ~/.bash_profile; FastTree -gtr -nt -gamma RR_DUP240.aln.nuc > RR_DUP240.aln.fasttree")
system("source ~/.bash_profile; FastTree -gtr -nt -gamma RR_DUP240_allAlleles.aln.nuc > RR_DUP240_allAlleles.aln.fasttree")
system("source ~/.bash_profile; FastTree -gamma RR_DUP240.aln.prot > tree_prot.fasttree")
tree_n <- read.newick("tree.fasttree")
tree_p <- read.newick("tree_prot.fasttree")
p <- ggtree(tree_p) + theme_tree2()  # standard topology

label_hypersensitive <- which(grepl("YPS1009|I14|YJM981|_RM", p$data$label)) # _ is important - otherwise get duplications (PRM8/9 has RM in it)
label_intermediate <- which(grepl("M22|BY|CLIB219|CLIB413", p$data$label))
label_resistant <- which(grepl("Y10|273164|YPS163|CBS2888|YJM978|YJM145|YJM454|PW5", p$data$label))
#label_all <- c(label_resistant, label_hypersensitive, label_medium); which(duplicated(label_all)) # test for duplicates - run once, just to make sure there are none
pheno_vector <- vector(length=nrow(p$data), mode="character")
pheno_vector[label_hypersensitive] <- "hypersensitive";
pheno_vector[label_intermediate] <- "intermediate";
pheno_vector[label_resistant] <- "resistant";
pheno_vector[which(pheno_vector == "")] <- NA
p$data$phenotype <- pheno_vector
p <- p + geom_tiplab(aes(color=phenotype), parse=TRUE) + scale_color_manual(values=c("#FF0000", "#FFCC33", "#00FF00"))
p
ggsave("tree_prot.pdf", width = 40, height = 20, units = "in")


#PRANK codon-aware alignment
# system("source ~/.bash_profile; prank -d=./RR_DUP240/RR_DUP240.fsa -o=RR_DUP240_r000001_ext001.codon.aln -codon -F -gaprate=0.00001 -gapext=0.01")
# system("source ~/.bash_profile; prank -d=./RR_DUP240/REF_DUP240.fsa -o=REF_DUP240_r0000001_ext0001.codon.aln -codon -F -gaprate=0.000001 -gapext=0.001")
# system("source ~/.bash_profile; prank -d=./RR_DUP240/REF_DUP240_subset.fsa -o=REF_DUP240_subset_r0000001_ext0001_no_F.codon.aln -codon -gaprate=0.000001 -gapext=0.001")

# Method 1 - doesn't work, but good to have as a backup way
# alignment_codon <- readDNAMultipleAlignment("REF_DUP240_subset_r0000001_ext0001_no_F.codon.aln.best.fas", "fasta")
# rownames(alignment_codon) <- paste(rownames(alignment_codon), "  ")
# write.phylip(alignment_codon, "codeml/RR_DUP240.codon.aln.nuc")

# Method 2
alignment_codon_ape <- ape::read.FASTA("REF_DUP240_subset_curated.fas")
names(alignment_codon_ape) <- paste(names(alignment_codon_ape), " ")
#names(alignment_codon_ape) <- paste0("s_", 1:18, "  ")
ape::write.dna(alignment_codon_ape, "codeml/RR_DUP240.codon.aln.nuc", nbcol=1, colsep=" \n", colw=1000000)

system("source ~/.bash_profile; FastTree -gtr -nt -gamma REF_DUP240_subset_curated.fas > codeml/tree.fasttree") # a tree is required for running PAML
# system("source ~/.bash_profile; cd codeml; codeml") # runs codeml from PAML


#### dN/dS Analysis by GARD Recombination Breakpoints ####
aln <- readDNAMultipleAlignment("REF_DUP240_subset_curated.fas", "fasta")
b1 <- c(320)
b2 <- c(89, 319)
b3 <- c(89, 319, 647)
b4 <- c(89, 266, 371, 647)
b5 <- c(62, 139, 266, 371, 646)
b6 <- c(62, 139, 265, 297, 378, 646)
b7 <- c(55, 108, 141, 266, 298, 378, 647)
bl <- list(b1, b2, b3, b4, b5, b6, b7)

# Round breakpoints to the nearest codon
bl <- lapply(bl, function(x) {return(round(x/3) * 3)}) # potentially change to floor() or ceiling()

# Append start/end positions
bl <- lapply(bl, function(x) {return(append(x, 0, after=0))})
bl <- lapply(bl, function(x) {return(append(x, 735))})
# check: lapply(bl, function(x) {return(x %% 3)})

ctl_body <- read_file("codeml/codeml.ctl")
ctl_list <- tstrsplit(ctl_body, split="\n")[4:33] #skip header

for (i in 1:7) {
  for (j in 1:(i+1)) {
    temp_aln <- DNAMultipleAlignment( aln, start=bl[[i]][j]+1, end=bl[[i]][j+1] )
    bin_aln <- as.DNAbin(temp_aln)
    filename <- paste0(i, "bp_sub", j)
    rownames(bin_aln) <- paste0(rownames(bin_aln), " ")
    ape::write.dna(bin_aln, paste0("codeml/split_aln/aln/", filename, ".nuc"), nbcol=1, colsep="  \n", colw=1000000)
    write.FASTA(bin_aln, paste0("codeml/split_aln/fasta/", filename, ".fsa"))
    system(paste0("source ~/.bash_profile; FastTree -gtr -nt -gamma codeml/split_aln/fasta/", filename, ".fsa > codeml/split_aln/tree/", filename, ".fasttree"))
    header <- paste0("      seqfile = ../aln/", filename, ".nuc\n      outfile = ../out/", filename, ".codeml\n     treefile = ../tree/", filename, ".fasttree\n")
    body <- append(ctl_list, tstrsplit(header, "\n"), after=0)
    write(paste0(body, collapse = "\n"), paste0("codeml/split_aln/ctl/", filename, ".ctl"))
  }
}

tree_1 <- read.newick("codeml/split_aln/tree/7bp_sub1.fasttree")
tree_2 <- read.newick("codeml/split_aln/tree/7bp_sub2.fasttree")
tree_3 <- read.newick("codeml/split_aln/tree/7bp_sub3.fasttree")
tree_4 <- read.newick("codeml/split_aln/tree/7bp_sub4.fasttree")
tree_5 <- read.newick("codeml/split_aln/tree/7bp_sub5.fasttree")
tree_6 <- read.newick("codeml/split_aln/tree/7bp_sub6.fasttree")
tree_7 <- read.newick("codeml/split_aln/tree/7bp_sub7.fasttree")
tree_8 <- read.newick("codeml/split_aln/tree/7bp_sub8.fasttree")

tree_list <- list(tree_1, tree_2, tree_3, tree_4, tree_5, tree_6, tree_7, tree_8)

i <- 1
for (tree in tree_list) {
  ggtree(tree) + theme_tree2() + geom_tiplab( aes(label=tstrsplit(label, "_")[[2]]) ) + geom_nodelab(aes(label=signif(100*as.numeric(label), 2)), nudge_x = .012)
  ggsave(paste0("~/Desktop/T27/tree_", i, ".pdf"), width = 7, height = 5.25)
  i <- i + 1
}

( p1 <- ggtree(tree_1) + 
    theme_tree2() + 
    geom_tiplab() + 
    #geom_text2(aes(subset = !isTip, label=signif(100*as.numeric(label)), 2)) +
    #geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 10)) +
    geom_nodelab(aes(label=signif(100*as.numeric(label), 2)), nudge_x = .01) + 
    coord_cartesian(xlim = c(0,1))
)

( p2 <- ggtree(tree_2) + theme_tree2() + geom_tiplab() )


# Fig. 3 - all KTD1 alleles visualization
RR_DUP240 <- readDNAStringSet(paste0(getwd(),"/DUP_YAR028W_allAlleles.fasta"))
#RR_DUP240 <- subseq(RR_DUP240, end=width(RR_DUP240)-3) # remove stop codon
names(RR_DUP240) <- paste(names(RR_DUP240), " ")# extra space required for conversion to Phylip below
nseq <- length(RR_DUP240) #number of DUP240s analyzed

# Multiple sequence alignment - may take a while
alignment_n <- DNAMultipleAlignment(msa(RR_DUP240, order="input"), use.names=TRUE) # convert to MultipleAlignment object
write.phylip(alignment_n, "RR_DUP240_allAlleles.aln.nuc")

# run "FastTree -gtr -nt -gamma RR_DUP240_allAlleles.aln.nuc > KTD1.fasttree"

# Visualize the tree
tree <- read.newick("KTD1.fasttree")
ggtree(tree) + theme_tree2() + geom_tiplab() #+ coord_cartesian(xlim = c(0, 0.06))

# meowmeow
RR_DUP240 <- readDNAStringSet(paste0(getwd(),"/DUP_YAR028W_allAlleles1.fasta"))
names(RR_DUP240) <- paste(names(RR_DUP240), " ")# extra space required for conversion to Phylip below
nseq <- length(RR_DUP240) #number of DUP240s analyzed

# Multiple sequence alignment - may take a while
alignment_n <- DNAMultipleAlignment(msa(RR_DUP240, order="input"), use.names=TRUE) # convert to MultipleAlignment object
writeXStringSet(as(unmasked(alignment_n), "XStringSet"), file="aln.fasta")


# T28 - YAR028W alleles gene tree
tree <- read.newick("RR_DUP240_allAlleles.aln.fasttree")
ggtree(tree) + theme_tree2() + geom_tiplab()# + coord_cartesian(xlim = c(0, 0.35))
#ggsave("YAR028W_alleles_geneTree.pdf", width = 9, height = 6)
ggsave("RR_DUP240_allAlleles.pdf", width = 40, height = 20, units = "in")

# Homoplasies removed
pchisq(2*(-1512.671767+1530.129935), df=2, lower.tail=FALSE)
