library(systemPipeR)

# sets working directory to the script location - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
script_dir <- getwd()
source("../utils.R")


#### INITIAL ROUND OF SEARCH: R0 aka REF, i.e. starting from just the reference DUPs and expanding the gene set ####
blast.dt <- fread("REFonly_vs_RR.out") #format: blastn ... -outfmt 6
colnames(blast.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.dt[, c("chr","strain") := tstrsplit(sseqid, split="_")]
blast.dt[sstart > send, c("sstart", "send") := .(send, sstart)]


blast.dt <- blast.dt[evalue < 0.1]
c.dt <- blast.dt[length > 650 & pident > 96, c("qseqid", "sstart", "send", "strain", "chr")]
c.dt[, c("h.strain", "h.chr") := .(strain, chr)]
setkey(c.dt, strain, chr, sstart, send)
overlaps <- foverlaps(blast.dt, c.dt, type="any")[, .(qseqid, h.strain, h.chr, length, sstart, send)]

#overlaps[is.na(qseqid), c("qseqid", "h.strain", "h.chr") := rep("NA", 107)]
#overlaps[, qseqid := paste(" Overlap Hit:", qseqid, h.strain, h.chr, length, sstart, send)]
#overlaps[, c("h.strain", "h.chr", "length", "sstart", "send") := rep(NULL, 4)]
#fwrite(overlaps, "overlaps.dat", col.names = FALSE, eol = "\n")

c.dt[, .N, by=qseqid]

# Identifying new/missed DFP's
overlaps <- foverlaps(blast.dt, c.dt, type="any")
seqbatch <- overlaps[length > 650 & is.na(h.chr), .(sseqid, i.sstart, i.send)][, length := i.send - i.sstart + 1]
setkey(seqbatch, length)
seqbatch <- seqbatch[, tail(.SD, 1), by = .(i.sstart)] # select for longest overlapping region - NOT a great method, but works in this case
seqbatch[, sseqid := paste0(sseqid, " -range ", i.sstart-100, "-", i.send+100)] #100bp flanks
fwrite(seqbatch[, .(sseqid)], "seqbatchREF.txt", col.names = FALSE, eol="\n")

# Melt c.dt to make a DUP table
o.melt <- c.dt[, .(qseqid, strain)]
o.melt <- dplyr::count(o.melt, qseqid, strain)
o.cast <- dcast(o.melt, strain ~ qseqid, value.var = "n", fill = 0)

# Reorder rows and columns
strain_index <- c(6,8,14,12,5,4,2,7,1,11,10,9,16,3,13,15)
o.cast[, index := strain_index]; setkey(o.cast, index); o.cast[, index := NULL]
newcolorder <- c("strain", "DFP1", "REF_UIP3", "REF_YAR028W", "REF_PRM9",
                 "REF_MST28", "REF_PRM8", "DFP3", "DFP4")
setcolorder(o.cast, newcolorder)


#### ROUND 1 (R1) ####
blast.dt <- fread("DFPandREF_R1_vs_RR.out") #format: blastn ... -outfmt 6
colnames(blast.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.dt[, c("chr","strain") := tstrsplit(sseqid, split="_")]
blast.dt[sstart > send, c("sstart", "send") := .(send, sstart)]
blast.dt[, "index" := 1:.N]
blast.dt <- blast.dt[evalue < 0.1]

c.dt <- blast.dt[length > 650 & pident > 96]
c.dt <- rbind(c.dt, blast.dt[index %in% c(47, 52, 57, 62, 440, 452)]) #check for qseqid == "DFP1" & length > 650 (and == "REF_PRM8" & 650)
c.dt <- c.dt[, c("qseqid", "sstart", "send", "strain", "chr")]
c.dt[, c("h.strain", "h.chr") := .(strain, chr)]
setkey(c.dt, strain, chr, sstart, send)
overlaps <- foverlaps(blast.dt, c.dt, type="any")[, .(qseqid, h.strain, h.chr, length, sstart, send)]

#overlaps[is.na(qseqid), c("qseqid", "h.strain", "h.chr") := rep("NA", overlaps[is.na(qseqid), .N])]
#overlaps[, qseqid := paste(" Overlap Hit:", qseqid, h.strain, h.chr, length, sstart, send)]
#overlaps[, c("h.strain", "h.chr", "length", "sstart", "send") := rep(NULL, 4)]
#fwrite(overlaps, "overlaps_R1.dat", col.names = FALSE, eol = "\n")

# Identifying new/missed DFP's
overlaps <- foverlaps(blast.dt, c.dt, type="any")
seqbatch <- overlaps[length > 650 & is.na(h.chr), .(sseqid, i.sstart, i.send)][, length := i.send - i.sstart]
setkey(seqbatch, length)
seqbatch <- seqbatch[, tail(.SD, 1), by = .(i.sstart)] # select for longest overlapping region - NOT a great method, but works in this case
seqbatch[, sseqid := paste0(sseqid, " -range ", i.sstart-100, "-", i.send+100)] #100bp flanks - maybe not necessary
fwrite(seqbatch[, .(sseqid)], "seqbatchDFPandREF_R1.txt", col.names = FALSE, eol="\n")



#### ROUND 2 (R2) ####
blast.dt <- fread("DFPandREF_R2_vs_RR.out") #format: blastn ... -outfmt 6
colnames(blast.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.dt[, c("chr","strain") := tstrsplit(sseqid, split="_")]
blast.dt[sstart > send, c("sstart", "send") := .(send, sstart)]
blast.dt[, "index" := 1:.N]
blast.dt <- blast.dt[evalue < 0.1]

c.dt <- blast.dt[length > 650 & pident > 96]
c.dt <- rbind(c.dt, blast.dt[index %in% c(47, 52, 57, 62, 440, 452, 3714)]) #check for qseqid == "DFP1" & length > 650 (and == "REF_PRM8" & 650)
c.dt <- c.dt[, c("qseqid", "sstart", "send", "strain", "chr")]
c.dt[, c("h.strain", "h.chr") := .(strain, chr)]
setkey(c.dt, strain, chr, sstart, send)
#overlaps <- foverlaps(blast.dt, c.dt, type="any")[, .(qseqid, h.strain, h.chr, length, sstart, send)]

#overlaps[is.na(qseqid), c("qseqid", "h.strain", "h.chr") := rep("NA", overlaps[is.na(qseqid), .N])]
#overlaps[, qseqid := paste(" Overlap Hit:", qseqid, h.strain, h.chr, length, sstart, send)]
#overlaps[, c("h.strain", "h.chr", "length", "sstart", "send") := rep(NULL, 4)]
#fwrite(overlaps, "overlaps_R1.dat", col.names = FALSE, eol = "\n")

# Identifying new/missed DFP's
overlaps <- foverlaps(blast.dt, c.dt, type="any")
# No hits > 650!
seqbatch <- overlaps[length > 500 & is.na(h.chr), .(sseqid, i.sstart, i.send)][, length := i.send - i.sstart]
setkey(seqbatch, length)
seqbatch <- seqbatch[, tail(.SD, 1), by = .(i.sstart)] # select for longest overlapping region - NOT a great method, but works in this case
seqbatch[, sseqid := paste0(sseqid, " -range ", i.sstart-500, "-", i.send+500)] #500 bp flanks
fwrite(seqbatch[, .(sseqid)], "seqbatchDFPandREF_R2.txt", col.names = FALSE, eol="\n")

# Any hits overlapping with DFP19 in CLIB219? - Yes, but no full-length ones! Need to extract ORFs
just19 <- data.table(qseqid = "DFP19", sstart = 175684, send = 176349, strain = "CLIB219", chr = "chrI", h.strain = "CLIB219", h.chr = "chrI")
setkey(just19, strain, chr, sstart, send)
overlaps19 = foverlaps(blast.dt, just19, type = "any")
overlaps19[!is.na(qseqid)]


# EXTRACTING ORFs AT R2

# 0: Pull sequences with ranges (hitL - 1000), (hitR + 1000) [safe 1kb flank]
# -- these are in overlaps[is.na(h.chr)]

seqbatch_na <- overlaps[is.na(h.chr), .(sseqid, i.sstart, i.send)][, length := i.send - i.sstart]
setkey(seqbatch_na, length, i.sstart)
seqbatch_na[, sseqid := paste0(sseqid, " -range ", i.sstart-1000, "-", i.send+1000)]
fwrite(seqbatch_na[, .(sseqid)], "seqbatch_na_DFPandREF_R2.txt", col.names = FALSE, eol="\n")

# 1: Extract hits with BLAST+, import a FASTA file of pulled hits (lots of hits, most are short)
seqs <- readDNAStringSet("hits_na_DFPandREF_R2.fasta")

# 2: Processing (remove duplicates, etc.)
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE ) #fixes names of seqs
seqs <- seqs[!duplicated(names(seqs))] # alternative: unique(seqs), but this doesn't seem right

# 3: Run ORF caller
predicted_ORFs <- unlist(systemPipeR::predORF(seqs, n='all', strand='both')) #alt: n=1 for the longest ORF
predicted_ORFs_p <- predicted_ORFs[which(width(predicted_ORFs) >= 600)]
names(predicted_ORFs_p) <- seqnames(predicted_ORFs_p)
predicted_ORFs_p <- predicted_ORFs[ order(predicted_ORFs, names(predicted_ORFs), strand(predicted_ORFs), start(predicted_ORFs), end(predicted_ORFs)) ]
predicted_ORFs_p <- subset(predicted_ORFs_p, end(predicted_ORFs_p) >= 1000 & start(predicted_ORFs_p) <= 1000) #excludes ORFs not crossing pos 1000

# 4: Recover sequences
ORFseqs <- getSeq(seqs, predicted_ORFs_p)
names(ORFseqs) <- paste0(names(predicted_ORFs_p),"_",start(predicted_ORFs_p),"..",end(predicted_ORFs_p),"(",strand(predicted_ORFs_p),")") # assigns an appropriate ORF name based on location on the trimmed node

# 5: Cross-reference by BLASTing them against reference DUPS
writeXStringSet(ORFseqs, "new_ORFseqs_R3.fasta")
blast.dt <- fread("new_ORFseqs_R3.out") #format: blastn ... -outfmt 6
colnames(blast.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

print(blast.dt[, head(.SD, 1), by = .(qseqid)][, head(.SD, 1), by = .(pident)])



#### ROUND 3 (R3) ####
blast.dt <- fread("DFPandREF_R3_vs_RR.out") #format: blastn ... -outfmt 6
colnames(blast.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast.dt[, c("chr","strain") := tstrsplit(sseqid, split="_")]
blast.dt[sstart > send, c("sstart", "send") := .(send, sstart)]
blast.dt[, "index" := 1:.N]
blast.dt <- blast.dt[evalue < 0.1]

c.dt <- blast.dt[length > 650 & pident > 96]
c.dt <- rbind(c.dt, blast.dt[index %in% c(47, 52, 57, 62, 440, 452, 3714, 6864, 6878)]) #check blast.dt for qseqid == "DFP1" & length > 650 (and == "REF_PRM8" & 650) (and == "DFP26" & 650)
# DFP19 is taken care of.
# it is very nice that BLAST is deterministic - gives same hits in the same order
c.dt <- c.dt[, c("qseqid", "sstart", "send", "strain", "chr")]
c.dt[, c("h.strain", "h.chr") := .(strain, chr)]
setkey(c.dt, strain, chr, sstart, send)

overlaps <- foverlaps(blast.dt, c.dt, type="any")
# No hits > 650! overlaps[is.na(h.chr) & length > 650] == empty dt
# overlaps[is.na(h.chr) & length > 500] - only 5 hits
#seqbatch <- overlaps[length > 500 & is.na(h.chr), .(sseqid, i.sstart, i.send)][, length := i.send - i.sstart]
#setkey(seqbatch, length)
#seqbatch <- seqbatch[, tail(.SD, 1), by = .(i.sstart)] # select for longest overlapping region - NOT a great method, but works in this case
#seqbatch[, sseqid := paste0(sseqid, " -range ", i.sstart-500, "-", i.send+500)] #500 bp flanks
#fwrite(seqbatch[, .(sseqid)], "seqbatchDFPandREF_R2.txt", col.names = FALSE, eol="\n")

# EXTRACTING ORFs AT R3

# 0: Pull sequences with ranges (hitL - 1000), (hitR + 1000) [safe 1kb flank]
# -- these are in overlaps[is.na(h.chr)]

seqbatch_na <- overlaps[is.na(h.chr), .(sseqid, i.sstart, i.send)][, length := i.send - i.sstart]
setkey(seqbatch_na, length, i.sstart)
seqbatch_na[, sseqid := paste0(sseqid, " -range ", i.sstart-1000, "-", i.send+1000)]
fwrite(seqbatch_na[, .(sseqid)], "seqbatch_na_DFPandREF_R3.xt", col.names = FALSE, eol="\n")

# 1: Extract hits with BLAST+, import a FASTA file of pulled hits (lots of hits, most are short)
seqs <- readDNAStringSet("hits_na_DFPandREF_R3.fasta")

# 2: Processing (remove duplicates, etc.)
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE ) #fixes names of seqs
seqs <- seqs[!duplicated(names(seqs))] 

# 3: Run ORF caller
predicted_ORFs <- unlist(systemPipeR::predORF(seqs, n='all', strand='both')) #alt: n=1 for the longest ORF
predicted_ORFs_p <- predicted_ORFs[which(width(predicted_ORFs) >= 600)]
names(predicted_ORFs_p) <- seqnames(predicted_ORFs_p)
# now we go back to predicted_ORFs
predicted_ORFs_p <- predicted_ORFs[ order(predicted_ORFs, names(predicted_ORFs), strand(predicted_ORFs), start(predicted_ORFs), end(predicted_ORFs)) ]
predicted_ORFs_p <- subset(predicted_ORFs_p, end(predicted_ORFs_p) >= 1000 & start(predicted_ORFs_p) <= 1000) #excludes ORFs not crossing pos 1000

# 4: Recover sequences
ORFseqs <- getSeq(seqs, predicted_ORFs_p)
names(ORFseqs) <- paste0(names(predicted_ORFs_p),"_",start(predicted_ORFs_p),"..",end(predicted_ORFs_p),"(",strand(predicted_ORFs_p),")") # assigns an appropriate ORF name based on location on the trimmed node

# 5: Cross-reference by BLASTing them against reference DUPS
writeXStringSet(ORFseqs, "new_ORFseqs_R3.fasta")
# run blastn via command console
blast.dt <- fread("new_ORFseqs_R3.out") #format: blastn ... -outfmt 6
colnames(blast.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Processing
blast.dt[, c("chr","strain","rcoords") := tstrsplit(qseqid, split="_")]
blast.dt[, c("strain", "region") := tstrsplit(strain, split=":")]
blast.dt[, c("rcoords", "strand") := tstrsplit(rcoords, split="(", fixed=TRUE)]
blast.dt[, "strand" := fifelse(substr(strand, 1, 1) == "+", "plus", "minus")]
blast.dt[, c("rstart", "rend") := tstrsplit(rcoords, split="..", fixed=TRUE)] 
blast.dt[, c("rstart", "rend") := lapply(.SD, as.numeric), .SDcols = c("rstart", "rend")]
blast.dt[, "offset" := tstrsplit(region, split="-")[1]]; blast.dt[, "offset" := as.numeric(offset) - 1]
blast.dt[, c("rstart", "rend") := lapply(.SD, function(x){return(x+offset)}), .SDcols = c("rstart", "rend")]
blast.dt[, c("sstart", "send") := .(rstart, rend)]
#setkey(blast.dt, sstart, send)

# Overlaps - looking for "signatures" of DUPs, excluding any hits to reference DUPs (in c.dt)
overlaps <- foverlaps(blast.dt, c.dt, type="any")
overlaps[, "ORFlength" := rend - rstart + 1]
overlaps[is.na(h.chr) & length > 200]
overlaps[is.na(h.chr) & length > 100, head(.SD, 1), by = .(pident)]
overlaps[is.na(h.chr) & length > 100, head(.SD, 1), by = .(rstart)] #not ideal, but likelihood of having same rstart in different strains is low
overlaps[is.na(h.chr) & length > 100, head(.SD, 1)[, head(.SD, 1), by = .(rstart)], by = .(strain)]

# Verifying that ORFs start/end positions are correct (by translation)
# !!THE CODE BELOW MUST BE RUN ON A blast.dt VERSION FROM THE BEGINNING OF R3 CODE SECTION!!
#y029_ind <- blast.dt[qseqid == "REF_YAR029W" & pident > 91, head(.SD, 1), by = .(sseqid)]$index
y029_ind <- blast.dt[qseqid == "REF_YAR029W" & length > 200]$index #same as above, but note two hits from YPS163
c.dt <- rbind(c.dt, blast.dt[index %in% y029_ind]) #add YAR029W hits, ad-hoc
c.dt <- c.dt[, c("qseqid", "sstart", "send", "strain", "chr")]
c.dt[, c("h.strain", "h.chr") := .(strain, chr)]
setkey(c.dt, strain, chr, sstart, send)

blast.dt[, "batch" := paste0(chr, '_', strain, ' -range ', rstart, '-', rend, ' -strand ', strand)]
fwrite(blast.dt[, .(batch)], "seqbatch_newORFseqs_R3.txt", col.names = FALSE, eol="\n")
ORFs_verify <- readDNAStringSet("newORFseqs_R3_verification.fasta")
table(width(ORFs_verify)/3 == width(translate(ORFs_verify))) # true

# Verifying c.dt contains ORFs
c.dt[, "batch" := paste0(chr, '_', h.strain, ' -range ', sstart, '-', send)]
# TO-Do

blast.dt[sstart > send, c("sstart", "send") := .(send, sstart)] # s -> q?
blast.dt[, "index" := 1:.N]
#blast.dt <- blast.dt[evalue < 0.1]

print(blast.dt[, head(.SD, 1), by = .(qseqid)])
print(blast.dt[, head(.SD, 1), by = .(qseqid)][, head(.SD, 1), by = .(pident)])
blast.dt[length > 650] #3 hits, all mapping to the same location
blast.dt[length > 600] #6 hits, same as above
blast.dt[length > 500] #62 hits
blast.dt[length > 500, print(.SD), by = .(qseqid)]

#yar029w_add <- data.table(qseqid = "REF_YAR029W", sstart = 175434, send = 176009, strain = "YPS163", chr = "chrI")
# c.dt <- rbind(c.dt, yar029w_add)


# DUP table
o.melt <- c.dt[, .(qseqid, strain)]
o.melt <- dplyr::count(o.melt, qseqid, strain)
o.cast <- dcast(o.melt, strain ~ qseqid, value.var = "n", fill = 0)

# Reorder rows and columns
strain_index <- c(6,8,14,12,5,4,2,7,1,11,10,9,16,3,13,15)
o.cast[, index := strain_index]; setkey(o.cast, index); o.cast[, index := NULL]
newcolorder <- c("strain", "DFP1", "REF_UIP3", "REF_YAR028W", "REF_YAR029W", "REF_PRM9",
                 "REF_MST28", "REF_PRM8", "DFP3", "DFP4", paste0("DFP", 11:26))
setcolorder(o.cast, newcolorder)
fwrite(o.cast, "DUP240_table-complete.csv") #export


# Translating c.dt to identify *'d DUPs 
c.dt[, "batch" := paste0(chr, '_', strain, ' -range ', sstart, '-', send)]
fwrite(c.dt[, .(batch)], "seqbatch_c-dt.txt", col.names = FALSE, eol="\n")
# extract with blastdbcmd here

seqs <- readDNAStringSet("c.dt.fasta")
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE) #fixes names of seqs
predicted_ORFs <- unlist(systemPipeR::predORF(seqs, n=1, strand='both')) #longest ORF, but for both strands

c.dt[, 'batch' := NULL]

predORF <- as.data.table(predicted_ORFs)[strand == "+"]
predORF[, 'orflengthP' := width]; predORF[, 'orflengthM' := as.data.table(predicted_ORFs)[strand == "-"]$width]
predORF[, c("start", "end", "width", "strand", "subject_id", "inframe2end") := rep(NULL, 6)]
predORF[, c("chr", "start") := tstrsplit(seqnames, split = ':')]
predORF[, c("chr", "strain") := tstrsplit(chr, split = '_')]
predORF[, c("start", "end") := tstrsplit(start, split = '-')]
predORF[, c("seqnames") := NULL]
predORF[, "max_orflength" := pmax(orflengthP, orflengthM)] #parallel max - nice!
setkey(predORF, chr, strain, start)

setkey(c.dt, chr, strain, sstart, send)
c.dt$sstart == predORF$start # check for same order
c.dt[, 'max_orflength' := predORF[, max_orflength]]

c.dt[qseqid != "REF_YAR029W" & max_orflength < 650]
c.dt[, "batch" := paste0(chr, '_', strain, ' -range ', sstart-1000, '-', send+1000)]
fwrite(c.dt[qseqid != "REF_YAR029W" & max_orflength < 650, .(batch)], "seqbatch_c-dt_disrupted.txt", col.names = FALSE, eol="\n")


# Translating c.dt to identify *'d DUPs - 400 BP FLANKS VERSION
c.dt[, "batch" := paste0(chr, '_', strain, ' -range ', sstart-400, '-', send+400)]
fwrite(c.dt[, .(batch)], "seqbatch_c-dt400.txt", col.names = FALSE, eol="\n")
# extract with blastdbcmd here

seqs <- readDNAStringSet("c.dt400.fasta")
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE) #fixes names of seqs
predicted_ORFs <- unlist(systemPipeR::predORF(seqs, n=1, strand='both')) #longest ORF, but for both strands

predORF <- as.data.table(predicted_ORFs)[strand == "+"]
predORF[, 'orflengthP' := width]; predORF[, 'orflengthM' := as.data.table(predicted_ORFs)[strand == "-"]$width]
predORF[, c("start", "end", "width", "strand", "subject_id", "inframe2end") := rep(NULL, 6)]
predORF[, c("chr", "start") := tstrsplit(seqnames, split = ':')]
predORF[, c("chr", "strain") := tstrsplit(chr, split = '_')]
predORF[, c("start", "end") := tstrsplit(start, split = '-')]
predORF[, "seqnames" := NULL]
predORF[, "max_orflength" := pmax(orflengthP, orflengthM)] #parallel max - nice!
setkey(predORF, chr, strain, start)

setkey(c.dt, chr, strain, sstart, send)
c.dt$sstart-400 == predORF$start # check for same order
c.dt[, 'max_orflength' := predORF[, max_orflength]]

c.dt[qseqid != "REF_YAR029W" & max_orflength < 700]
fwrite(c.dt[qseqid != "REF_YAR029W" & max_orflength < 700, .(batch)], "seqbatch_c-dt_disrupted400.txt", col.names = FALSE, eol="\n")

# check the set of lengths
c2.dt <- c.dt[, head(.SD, 1), by = .(qseqid, max_orflength)]
setkey(c2.dt, qseqid, max_orflength)


#### Known arrays orthogonal cross-check ####
#array.dt <- fread("chrI_CDC15-YAT1.out")
array.dt <- fread("chrVII_ERV14-TYW3.out")

colnames(array.dt) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
setkey(array.dt, qseqid, sseqid)
array.dt <- array.dt[order(qseqid, sseqid)]
array.dt[1:16, sseqid] == array.dt[16 + 1:16, sseqid] #check for same order
array.dt[1:16, leftend := pmin(sstart, send)]
array.dt[1:16, rightend := array.dt[17:32, pmax(sstart, send)]]
array.dt <- array.dt[1:16][, qseqid := NULL]
array.dt[, "batch" := paste0(sseqid, ' -range ', leftend, '-', rightend)]
# fwrite(array.dt[, .(batch)], "seqbatch_CDC15-YAT1_flanks.txt", col.names = FALSE, eol="\n")
fwrite(array.dt[, .(batch)], "seqbatch_ERV14-TYW3_flanks.txt", col.names = FALSE, eol="\n")

# blastdbcmd-extract here
# seqs <- readDNAStringSet("RR_chrI_arrays.fasta")
seqs <- readDNAStringSet("RR_chrVII_arrays.fasta")
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE) #fixes names of seqs
predicted_ORFs <- unlist(systemPipeR::predORF(seqs, n='all', strand='both'))

predORF <- as.data.table(predicted_ORFs)
predORF[, "roffset" := tstrsplit(seqnames, split=":", fixed = TRUE)[[2]]]
predORF[, "roffset" := tstrsplit(roffset, split="-")[[1]]]
predORF[, "roffset" := as.numeric(roffset) - 1]
predORF[, "chr" := tstrsplit(seqnames, split = ':')[[1]] ]
predORF[, c("chr", "strain") := tstrsplit(chr, split = '_')]
predORF[, "start" := start + roffset]
predORF[, "end" := end + roffset]
predORF[strand == "-", c("start", "end") := .(end, start)] 
setkey(predORF, start, end)

predORF <- predORF[width > 650] #length >650 subset
predORF <- rbind(predORF[strand == "+", head(.SD, 1), by = end], predORF[strand == "-", tail(.SD, 1), by = end])
predORF[, c("seqnames", "strand", "subject_id", "inframe2end", "roffset") := NULL]
predORF[start > end, c("start", "end") := .(end, start)] 
setcolorder(predORF, c('chr', 'strain', 'width', 'start', 'end'))

c.dt <- fread("c.dt"); colnames(c.dt) <- c("qseqid", "start", "end", "strain", "chr")
setkey(c.dt, strain, chr, start, end)
overlaps <- foverlaps(predORF, c.dt, type="any")

#overlaps.chrI <- copy(overlaps)
#overlaps.chrVII <- copy(overlaps)

overlaps.chrI[is.na(qseqid)]
overlaps.chrVII[is.na(qseqid) & width != 822]

fwrite(overlaps.chrI, "chrI_array_ORFs_650bpCutoff.csv")
fwrite(overlaps.chrVII, "chrVII_array_ORFs_650bpCutoff.csv")


#### EXACT ORFs COORDINATES IN c.dt ####
c.dt <- fread("c.dt"); colnames(c.dt) <- c("qseqid", "sstart", "send", "strain", "chr")
c.dt[, "swidth" := send - sstart + 1]
c.dt[, "batch" := paste0(chr, '_', strain, ' -range ', sstart-400, '-', send+400)]
fwrite(c.dt[, .(batch)], "seqbatch_c-dt400.txt", col.names = FALSE, eol="\n")
# extract with blastdbcmd here

seqs <- readDNAStringSet("c.dt400.fasta")
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE) #fixes names of seqs
predicted_ORFs <- unlist(systemPipeR::predORF(seqs, n=1, strand='both')) #longest ORF, but for both strands

predORF <- as.data.table(predicted_ORFs)[strand == "+"]
predORF[, c("startP", "endP", "orflengthP") := .(start, end, width)]
predORF[, c("startM", "endM", "orflengthM") := as.data.table(predicted_ORFs)[strand == "-", .(start, end, width)]]
predORF[, c("start", "end", "width", "strand", "subject_id", "inframe2end") := rep(NULL, 6)]
predORF[, c("chr", "start") := tstrsplit(seqnames, split = ':')]
predORF[, c("chr", "strain") := tstrsplit(chr, split = '_')]
predORF[, c("start", "end") := tstrsplit(start, split = '-')]
predORF[, "seqnames" := NULL]
predORF[, "max_orflength" := pmax(orflengthP, orflengthM)] #parallel max - nice!
predORF[, 'strand' := ifelse(max_orflength == orflengthP, '+', '-')]
setkey(predORF, chr, strain, start, end)


setkey(c.dt, chr, strain, sstart, send)
c.dt$sstart-400 == predORF$start # check for same order
# c.dt[, index := 1:.N]
# ind029 <- c.dt[qseqid == "REF_YAR029W", index]
# predORF[ind029, c("max_orflength", "strand") := .(orflengthP, '+')]

c.dt[, c('max_orflength', 'strand') := predORF[, .(max_orflength, strand)]]
c.dt[strand == '-', c("ostart", "oend") := predORF[strand == '-', .(startM, endM)]]
c.dt[strand == '+', c("ostart", "oend") := predORF[strand == '+', .(startP, endP)]]
c.dt[, "roffset" := sstart - 401]
c.dt[, c("ostart", "oend") := .(ostart + roffset, oend + roffset)]

c.dt[, "orfbatch" := paste0(chr, '_', strain, ' -range ', ostart, '-', oend, ' -strand ', fifelse(strand == '-', 'minus', 'plus'))]
fwrite(c.dt[, .(orfbatch)], "seqbatch_c-dt400orfs.txt", col.names = FALSE, eol="\n")
#extract here

orfs <- as.data.table(readDNAStringSet("c.dt400orfs.fasta"))
c.dt[, "seq" := unlist(orfs[, x])]
c.dt[, c("batch", "orfbatch", "roffset") := rep(NULL, 3)]
c.dt[, "intact" := fifelse(swidth > max_orflength, 'no', 'yes')]
setcolorder(c.dt, c("qseqid", "strain", "chr", "strand", "sstart", "send", "swidth", "intact", "ostart", "oend", "max_orflength", "seq"))
c.dt[intact == "no"]

setkey(c.dt, qseqid, strain, chr, sstart)
fwrite(c.dt, "c.dt.csv")


# #YAR029W gets special treatment because there's a slightly longer ORF on the both strands
# ### TO-DO ##
# predicted_ORFs_all <- unlist(systemPipeR::predORF(seqs, n='all', strand='both')) #longest ORF, but for both strands
# predORF <- as.data.table(predicted_ORFs_all)[strand == "+"]
# predORF[, c("startP", "endP", "orflengthP") := .(start, end, width)]
# predORF[, c("chr", "start") := tstrsplit(seqnames, split = ':')]
# predORF[, c("chr", "strain") := tstrsplit(chr, split = '_')]
# predORF[, c("start", "end") := tstrsplit(start, split = '-')]
# predORF[, c("start", "end") := .(as.numeric(start), as.numeric(end))]
# predORF[, c("startP", "endP") := .(startP + start - 401, endP + start - 401)]
# #next up: find overlaps w/ 029W "start" and "end"
# 
# setkey(c.dt, strain, chr, sstart, send)
# setkey(predORF, strain, chr, startP, endP)
# overlaps <- foverlaps(predORF, c.dt[qseqid == "REF_YAR029W"], type="any")
# overlaps[order(-orflengthP)][!is.na(qseqid), head(.SD, 1), by = .(ostart)] #gives weird long orfs



#### ATTEMPT #2 AT YAR029W ####
c.dt[qseqid == "REF_YAR029W"]
c.dt[, "orfbatch" := paste0(chr, '_', strain, ' -range ', sstart, '-', sstart+600)] #600 bp flank - covers 576 bp max orf, expect everything else to be shorter
fwrite(c.dt[qseqid == "REF_YAR029W", .(orfbatch)], "seqbatch_029W_600bp.txt", col.names = FALSE, eol="\n")
#extract with blastdbcmd

seqs <- readDNAStringSet("c.dt029W.fasta")
names(seqs) <- sapply( names(seqs), function(x) { return(substr(x, 1, nchar(x)-1)) }, USE.NAMES=FALSE) #fixes names of seqs
predicted_ORFs_all <- unlist(systemPipeR::predORF(seqs, n='all', strand='sense')) #plus strand
predORF <- as.data.table(predicted_ORFs_all)
predORF[, c("startP", "endP", "orflengthP") := .(start, end, width)]
predORF[, c("chr", "start") := tstrsplit(seqnames, split = ':')]
predORF[, c("chr", "strain") := tstrsplit(chr, split = '_')]
predORF[, c("start", "end") := tstrsplit(start, split = '-')]
predORF[, c("start", "end") := .(as.numeric(start), as.numeric(end))]
# predORF[, c("startP", "endP") := .(startP + start - 1, endP + start - 1)]


predORF[startP == 1] #I14 verified CAA -> TAA (nonsense)

#check same order
predORF[, c("startP", "endP") := .(startP + start - 1, endP + start - 1)]
setkey(predORF, strain, chr, startP, endP)
c.dt[qseqid == "REF_YAR029W", sstart] == predORF[startP == start, start] # prior to reassignment of startP/endP, the condition is 'startP == 1'
c.dt[qseqid == "REF_YAR029W", strand := '+']
c.dt[qseqid == "REF_YAR029W", c("ostart", "oend", "max_orflength") := predORF[startP == start, .(startP, endP, orflengthP)]]

c.dt[, "orfbatch" := paste0(chr, '_', strain, ' -range ', ostart, '-', oend, ' -strand ', fifelse(strand == '-', 'minus', 'plus'))]
fwrite(c.dt[qseqid == "REF_YAR029W", .(orfbatch)], "seqbatch_029Worfs.txt", col.names = FALSE, eol="\n")
#extract

orfs <- as.data.table(readDNAStringSet("c.dt029Worfs.fasta"))
c.dt[qseqid == "REF_YAR029W", seq := unlist(orfs[, x])]
c.dt[, orfbatch := NULL]
c.dt[qseqid == "REF_YAR029W"]
c.dt[qseqid == "REF_YAR029W", max_orflength] == c.dt[qseqid == "REF_YAR029W", nchar(seq)] #check

c.dt[, "intact" := fifelse(swidth > max_orflength, 'no', 'yes')]
setcolorder(c.dt, c("qseqid", "strain", "chr", "strand", "sstart", "send", "swidth", "intact", "ostart", "oend", "max_orflength", "seq"))
c.dt[intact == "no"]

setkey(c.dt, qseqid, strain, chr, sstart)
fwrite(c.dt, "c.dt.new.csv")
