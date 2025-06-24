## dada2 QAQC of Oorc fecal fastq sequences
## 1/3/2022 - modified 5/27/25 by mball for WADE-003 Arctic Predator Diets
## Amy Van Cise

# Set library path to /gscratch (replace with your full path)
.libPaths("/gscratch/coenv/mball/WADE003-arctic-pred/R_libs")

# Create the directory if it doesn't exist
dir.create("/gscratch/coenv/mball/WADE003-arctic-pred/R_libs", recursive = TRUE, showWarnings = FALSE)

# Verify the path is set correctly
.libPaths()

### set up working environment

devtools::install_github("benjjneb/dada2", ref="v1.16", lib = .libPaths()[1])
BiocManager::install("dada2", lib = .libPaths()[1])
BiocManager::install("S4Vectors")

library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)

diet.seqs.file <- "/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata/16SP2/"
taxref <- "/gscratch/coenv/mball3/WADE003-arctic-pred/16S_Arctic_predator_reference_database_05_2025.fasta"


### read fastq files in working directory
fnFs <- sort(list.files(diet.seqs.file, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(diet.seqs.file, pattern="_R2_001.fastq", full.names = TRUE))
sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")

### vizualize read quality profiles
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

### Name filtered files in filtered/subdirectory
filtFs <- file.path(diet.seqs.file, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(diet.seqs.file, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


### Filter and Trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
### Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### Learn Error Rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

#plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

### Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

### Construct sequence table
seqtab <- makeSequenceTable(mergers)

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names


### Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, taxref, tryRC = TRUE, minBoot = 95, outputBootstraps = FALSE)

### Save data
save(seqtab.nochim, freq.nochim, track, taxa, file = "WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
