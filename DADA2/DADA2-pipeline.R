# DADA2 Pipeline

# ------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") 

install.packages(DECIPHER)
# ------------------------------------------------------------------
# Load in the Library; Sets path
# ------------------------------------------------------------------
library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")

# identifies the path
path <- "../DADA2/MiSeq_SOP" 
list.files(path)

# ------------------------------------------------------------------
# Read in the Names of the fastq files; Match lists of F and R files
# ------------------------------------------------------------------
## file name format may differ

# Forward and reverse fastq filenames have format: 
# SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# ------------------------------------------------------------------
# Inspect Read QUality Profiles
# ------------------------------------------------------------------
# greyscale heatmap - frequency of each quality score at each base position
# green line - mean quality score at each position
# orange lines - quartiles of the quality score distribution 
# red line - scaled proportion of reads that extend to at least that posiotion (not super useful for illumina)

# Forward Reads - look good; could trim last nucleotides for quality control
plotQualityProfile(fnFs[1:2])

# Reverse reads -  significantly worse quality (esp. at the end; this is common for illumina)
## Though teh DADA2 algorithm incorporates quality information into the error model, its good to trim toimprove sensitivity
plotQualityProfile(fnRs[1:2])

# ------------------------------------------------------------------
# Filter and Trim
# ------------------------------------------------------------------
# To determine where to trim, use the x axis and the quality of the reads
## READS STILL MUST OVERLAP AFTER TRUNCATION

# Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtering; standard filtering parameters used here
## maxEE sets the maximum number of 'expected errors' allowed in a read (better than averaging quality scores)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160)#truncate low quality tails, 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


# ------------------------------------------------------------------
# Learn the Error Rates
# ------------------------------------------------------------------
# parametric error model (err); learnErrors literally learns errors from the data
## begins by guessing the max possible error rates in the data, then narrows until error and sample composition converge (???)
### The maximum errors for a data set = error rates if only the most abundant sequence is correct and all the rest are errors

# Forward Errors 
errF <- learnErrors(filtFs, multithread=TRUE)
## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

# Reverse Errors
errR <- learnErrors(filtRs, multithread=TRUE)
## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

# It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
  plotErrors(errF, nominalQ=TRUE)
  ## shows the error rates for each possible transition (A→C, A→G, …)
  ## points - observed error rates for each consensus quality score (??)
  ## black line - estimated error rates after convergence of the machine-learning algorithm
  ## red line - error rates expected under the nominal definition of the Q-score
  ### blakc line should fit the points; red line should decrease with increased quality scores
  ### this is an example of everything looking reasonable; could proceed with confidenc it my data looks liek this
  
# ------------------------------------------------------------------
# Sample Inference
# ------------------------------------------------------------------
# Applies the core sample inference algorithm
  
# Forward Reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
  
# Reverse Reads
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Inspects the dada-class object created above
dadaFs[[1]]
  ## here, DADA2 algorithm inferred 128 'true' sequence variants from the 1979 unique sequences in the first sample
  ## help("dada-class") for more options 

# ------------------------------------------------------------------
# Merge Paired Reads
# ------------------------------------------------------------------
# Merges Forward and Reverse Reads - gives you the full denoised sequences
## By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, 
## and are identical to each other in the overlap region (but these conditions can be changed via function arguments)
###  Paired reads that do not exactly overlap are removed by mergePairs

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# ------------------------------------------------------------------
# Contruct Sequence Table
# ------------------------------------------------------------------
# Amplicon sequence variant table (ASV) table\

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## The sequence table is a matrix 
## with rows corresponding to (and named by) the samples, 
## columns corresponding to (and named by) the sequence variants. 
## This table contains 293 ASVs, and the lengths of our merged sequences
### Sequences that are much longer or shorter than expected may be the result of non-specific priming.
### remove by (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256])

# ------------------------------------------------------------------
# Remove Chimeras
# ------------------------------------------------------------------
# The core dada method corrects substitution and indel errors, but chimeras remain.
## Chimeric sequences are identified if they can be exactly reconstructed 
## by combining a left-segment and a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # Chimeras only account for ~4% of the merged sequence reads (when accounting for abundance of variants)

### If most of your reads were removed as chimeric, upstream processing may need to be revisited. 
### In almost all cases this is caused by primer sequences with ambiguous nucleotides 
### that were not removed prior to beginning the DADA2 pipeline.

# ------------------------------------------------------------------
# Reack Reads Through the Pipeline
# ------------------------------------------------------------------
# THis is a final check: tells you how many reads made it through each step

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) # may cause concern if you lose a LOT of reads at a single step, throughout


# ------------------------------------------------------------------
# Assign Taxonomy
# ------------------------------------------------------------------
# Uses a native implementation of the naive Bayesian classifier method 
## The assignTaxonomy function takes as input a set of sequences to be classified 
## and a training set of reference sequences with known taxonomy, 
## and outputs taxonomic assignments with at least minBoot bootstrap confidence

#### FASTAS (??) -- need additional file {silva_species_assignment_v132.fa.gz}
 taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Extension (??) -- need above file
 taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz") 

# Inspects taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
