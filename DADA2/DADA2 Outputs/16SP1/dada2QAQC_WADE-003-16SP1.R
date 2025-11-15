###dada2 QAQC of DF&G dDNA sequences
###AVC
###Fall 2024

# Set up environment -----------------------------------------------------------

library(dada2); packageVersion("dada2")
library(seqinr)
library(taxa)

path <- "M:/Arctic predator diet sequences/16SP1"

# Get sample names -------------------------------------------------------------

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))[-1]
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))[-1]

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnRs[7:8])

# Filter and trim --------------------------------------------------------------
#new subdirectory 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 40, truncLen=c(280,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

filtFnew <- sort(list.files(paste0(path, "/filtered"), pattern="_F_filt.fastq", full.names = TRUE))
filtRnew <- sort(list.files(paste0(path, "/filtered"), pattern="_R_filt.fastq", full.names = TRUE))
# learn error rate -------------------------------------------------------------

errF <- learnErrors(filtFnew, multithread=TRUE)
errR <- learnErrors(filtRnew, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)

# Infer samples ----------------------------------------------------------------

dadaFs <- dada(filtFnew, err=errF, multithread=TRUE)
dadaRs <- dada(filtRnew, err=errR, multithread=TRUE)

# Merge paired reads -----------------------------------------------------------

mergers <- mergePairs(dadaFs, filtFnew, dadaRs, filtRnew, verbose=TRUE)

# sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras --------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through pipeline -------------------------------------------------

getN <- function(x) sum(getUniques(x))
outnew <- out %>% as.data.frame() %>% rename("input" = 1) %>% filter(input > 10)
track <- cbind(outnew, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
dim(track)

# Assign taxonomy --------------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, "G:/My Drive/06 ADFG diet study/Data analysis/QAQC pipeline/16S_Arctic_predator_reference_database_2024.fasta", 
                       minBoot = 95, multithread=FALSE) 

# new method? 
# library(DECIPHER); packageVersion("DECIPHER")
# 
# dna <- DNAStringSet(getSequences(seqtab.nochim))
# trainingSet <- readDNAStringSet("G:/My Drive/06 ADFG diet study/Data analysis/QAQC pipeline/16S_salmon_groundfish_reference_database_2022.fasta")
# trainingSet <- as_taxon(trainingSet)
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors

# Save data --------------------------------------------------------------------

save(seqtab.nochim, taxa, file = "ADFG_16SP1_dada2output_95boot_102424.Rdata")
