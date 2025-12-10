## dada2 QAQC of Arctic Predator fastq sequences
## 1/3/2022 - modified 5/27/25 by mball for WADE-003 Arctic Predator Diets
## Amy Van Cise

# Set library path to /gscratch (replace with your full path)
.libPaths("/gscratch/coenv/mball/WADE003-arctic-pred/R_libs")

# Create the directory if it doesn't exist
dir.create("/gscratch/coenv/mball/WADE003-arctic-pred/R_libs", recursive = TRUE, showWarnings = FALSE)

# Verify the path is set correctly
.libPaths()

### set up working environment

# devtools::install_github("benjjneb/dada2", ref="v1.16", lib = .libPaths()[1])
#BiocManager::install("dada2", lib = .libPaths()[1], force = TRUE)
#BiocManager::install("S4Vectors")

library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)
library(dplyr)

diet.seqs.file <- "/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata/12SP1"
taxref <- "/gscratch/coenv/mball3/WADE003-arctic-pred/MURI_MFU_07_2025.fasta"
speciesref <- "/gscratch/coenv/mball3/WADE003-arctic-pred/12S-AddSpecies_11-25.fasta"

### read fastq files in working directory
fnFs <- sort(list.files(diet.seqs.file, pattern="_R1_001.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(diet.seqs.file, pattern="_R2_001.fastq(\\.gz)?$", full.names = TRUE))
sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")

### vizualize read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

### Name filtered files in filtered/subdirectory
filtFs <- file.path(diet.seqs.file, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(diet.seqs.file, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### Filter and Trim
## for 12S trimLeft = 30, truncLen=c(130, 130)
## for 16S trimLeft = 40, truncLen=c(280, 180)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 30, truncLen=c(130, 130),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose=TRUE)
# 
# # Check which files exist for forwards and reverse
# exists_F <- file.exists(filtFs)
# exists_R <- file.exists(filtRs)
# 
# # Optionally display missing files:
# cat("Missing forward filtered files:\n")
# print(filtFs[!exists_F])
# 
# cat("Missing reverse filtered files:\n")
# print(filtRs[!exists_R])
# 
# # Subset to only those samples where both filtered files exist
# keep <- exists_F & exists_R
# filtFs <- filtFs[keep]
# filtRs <- filtRs[keep]
# sample.names <- sample.names[keep]

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

plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

### Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 20, verbose=TRUE)

### Construct sequence table 
seqtab <- makeSequenceTable(mergers)

# # # FOR 12S!! ONLY
obect <- as.data.frame(colnames(seqtab)) %>%
  rename("sequence" = 1) %>%
  mutate(seq_length = nchar(sequence)) %>%
  mutate(new_sequence = substr(sequence, 1, 185)) %>%
  mutate(seq2_length = nchar(new_sequence))

colnames(seqtab) <- obect$new_sequence

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names


### Assign Taxonomy
taxam <- assignTaxonomy(seqtab.nochim, taxref, tryRC = TRUE, minBoot = 95)

getwd()
# Loads in Assign Taxonomy (taxa) table
#rdata <- load("WADE003-arcticpred_dada2_QAQC_12SP1_output-addSpecies-130;30-2.Rdata")

# Assign Species
genus.species <- assignSpecies(seqtab.nochim, speciesref)

# Ensures all taxonomic levels agree for allk rows within an species assignment
taxa <- as.data.frame(taxam) %>% # saves taxa as a dataframe
  rownames_to_column("ASV") %>%  # makes ASV rownames to a column
  filter(is.na(Species)) %>% # keeps only those in Species column that are NA
  left_join(as.data.frame(genus.species) %>% # adds genus.species df to tax_table by ASV column
              rownames_to_column("ASV"),  by = "ASV") %>%  
  unite(col = addSpecies, Genus.y,Species.y, sep = " ") %>%  # combines genus.species df columns into one and renames column addspecies
  ungroup() %>%  # ungroups by ASV
  mutate(addSpecies = case_when(addSpecies == "NA NA"~NA, TRUE~addSpecies)) %>% # changes addSpecies column to R NA if NA NA characters and to what is in addSpecies column when there is a species
  mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%  # changes NA to spp. in addSpecies column
  mutate(.grp = ifelse(is.na(addSpecies),  paste0("NA_grp_", row_number()),  addSpecies)) %>%  # creates new column where NAs are specified by rowname (NA_grp_5)
  group_by(.grp) %>%  # groups by column .grp
  mutate(Class = if (length(unique(Class)) > 1) NA else Class) %>% # makes sure column and .grp agree; if not -> NA
  mutate(Order = if (length(unique(Order)) > 1) NA else Order) %>% 
  mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>% 
  mutate(Genus.x = if (length(unique(Genus.x)) > 1) NA else Genus.x) %>% 
  ungroup() %>%  # ungroups
  dplyr::rename("Genus" = Genus.x, "Species" = addSpecies) %>%  # renames weird column names to names that make sense
  select(-Species.x) %>%  # removes Species.x column
  bind_rows(as.data.frame(taxa) %>% # makes taxa a df and adds its columns
              rownames_to_column("ASV") %>%  # adds the ASVs as a column instead of rownames
              filter(!is.na(Species))) %>%  # adds all rows back in 
  mutate(.grp = ifelse(is.na(Species),  # changes .grp column to ???
                       paste0("NA_grp_", row_number()),  Species)) %>%  
  group_by(.grp) %>%  # groups by .grp
  fill(Order, Family, Genus, .direction = "updown") %>% # fills in NAs with the previous or next non NA value (adds correct value where disagreements were)
  ungroup() %>%  # ungroups
  select(-.grp) %>%  # removes .grp column
  column_to_rownames("ASV") %>%  # puts the columns ASV back to rownames
  as.matrix() # transforms df back to matrix

### Save data
save(seqtab.nochim, freq.nochim, track, taxa, file = "WADE003-arcticpred_dada2_QAQC_12SP1_output.Rdata")