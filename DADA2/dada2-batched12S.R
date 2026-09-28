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

# devtools::install_github("benjjneb/dada2", force = TRUE, lib = .libPaths()[1])
# BiocManager::install("dada2", version = '3.20', lib = .libPaths()[1], force = TRUE)
# BiocManager::install("S4Vectors", version = '3.20')
# install.packages("seqinr") 

library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)
library(dplyr)
library(rentrez)
library(progress)
library(Biostrings)
library(xml2)

packageVersion("dada2")

diet.seqs.file.P1 <- "/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata/12SP1"
diet.seqs.file.P2 <- "/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata/12SP2"

taxref <- "/gscratch/coenv/mball3/WADE003-arctic-pred/MURI_MFU_07_2025.fasta"
speciesref <- "/gscratch/coenv/mball3/WADE003-arctic-pred/ADDSPECIES_MURI_MFU_07_2025.fasta"

### read fastq files in working directory
fnFs.P1 <- sort(list.files(diet.seqs.file.P1, pattern="_R1_001.fastq(\\.gz)?$", full.names = TRUE))
fnRs.P1 <- sort(list.files(diet.seqs.file.P1, pattern="_R2_001.fastq(\\.gz)?$", full.names = TRUE))
sample.names1.P1 <- sapply(strsplit(basename(fnFs.P1), "_"), `[`, 1)
sample.names2.P1 <- sapply(strsplit(basename(fnFs.P1), "_"), `[`, 2)
sample.names.P1 <- paste(sample.names1.P1, sample.names2.P1, sep = "_")

fnFs.P2 <- sort(list.files(diet.seqs.file.P2, pattern="_R1_001.fastq(\\.gz)?$", full.names = TRUE))
fnRs.P2 <- sort(list.files(diet.seqs.file.P2, pattern="_R2_001.fastq(\\.gz)?$", full.names = TRUE))
sample.names1.P2 <- sapply(strsplit(basename(fnFs.P2), "_"), `[`, 1)
sample.names2.P2 <- sapply(strsplit(basename(fnFs.P2), "_"), `[`, 2)
sample.names.P2 <- paste(sample.names1.P2, sample.names2.P2, sep = "_")



### vizualize read quality profiles
plotQualityProfile(fnFs.P1[1:2])
plotQualityProfile(fnRs.P1[1:2])

plotQualityProfile(fnFs.P2[1:2])
plotQualityProfile(fnRs.P2[1:2])


### Name filtered files in filtered/subdirectory
filtFs.P1 <- file.path(diet.seqs.file.P1, "filtered.P1", paste0(sample.names.P1, "_F_filt.fastq.gz"))
filtRs.P1 <- file.path(diet.seqs.file.P1, "filtered.P1", paste0(sample.names.P1, "_R_filt.fastq.gz"))
names(filtFs.P1) <- sample.names.P1
names(filtRs.P1) <- sample.names.P1

filtFs.P2 <- file.path(diet.seqs.file.P2, "filtered.P2", paste0(sample.names.P2, "_F_filt.fastq.gz"))
filtRs.P2 <- file.path(diet.seqs.file.P2, "filtered.P2", paste0(sample.names.P2, "_R_filt.fastq.gz"))
names(filtFs.P2) <- sample.names.P2
names(filtRs.P2) <- sample.names.P2

# Checks issue of vector length not equal to names attribute 
# length(fnFs)
# length(filtFs)
# length(sample.names)
# 
# # Check first few files have matching read counts
# # Count reads in first 5 forward files (divide lines by 4)
# for(i in 1:min(5, length(fnFs))) {
#   lines <- system(paste("zcat", fnFs[i], "| wc -l"), intern=TRUE)
#   reads <- as.numeric(lines) / 4
#   cat(basename(fnFs[i]), reads, "\n")
# }
# 
# for(i in 1:min(5, length(fnRs))) {
#   lines <- system(paste("zcat", fnRs[i], "| wc -l"), intern=TRUE)
#   reads <- as.numeric(lines) / 4
#   cat(basename(fnRs[i]), reads, "\n")
# }


### Filter and Trim
## for 12S trimLeft = 30, truncLen=c(130, 130)
out.P1 <- filterAndTrim(fnFs.P1, filtFs.P1, fnRs.P1, filtRs.P1, trimLeft = 30, truncLen=c(130, 130),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE, verbose=TRUE)

out.P2 <- filterAndTrim(fnFs.P2, filtFs.P2, fnRs.P2, filtRs.P2, trimLeft = 30, truncLen=c(130, 130),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE, verbose=TRUE)

### Dereplicate
derepFs.P1 <- derepFastq(filtFs.P1, verbose=TRUE)
derepRs.P1 <- derepFastq(filtRs.P1, verbose=TRUE)

derepFs.P2 <- derepFastq(filtFs.P2, verbose=TRUE)
derepRs.P2 <- derepFastq(filtRs.P2, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs.P1) <- sample.names.P1
names(derepRs.P1) <- sample.names.P1

names(derepFs.P2) <- sample.names.P2
names(derepRs.P2) <- sample.names.P2

### Learn Error Rates

# use for MiSeqs that do not bin error rates

dadaFs.lrn.P1 <- dada(derepFs.P1, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF.P1 <- dadaFs.lrn.P1[[1]]$err_out
dadaRs.lrn.P1 <- dada(derepRs.P1, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR.P1 <- dadaRs.lrn.P1[[1]]$err_out

plotErrors(dadaFs.lrn.P1[[1]], nominalQ=TRUE)

# use for MiSeqs that do bin error rates

# Checks the current dada2 package version to use makeBinnedQualErrfun()
packageVersion("dada2") # package version should be v1.34.0 or newer

### Learn Error Rates (BINNED QUALITY SCORES METHOD using pre-learned approach)
novaBinnedErrfun <- makeBinnedQualErrfun(c(2, 12, 24, 38))

# Learn error rates FIRST from filtered reads (using binned error function) 
errF.P2 <- learnErrors(filtFs.P2, errorEstimationFunction = novaBinnedErrfun, multithread=TRUE, verbose = TRUE)

errR.P2 <- learnErrors(filtRs.P2, errorEstimationFunction = novaBinnedErrfun, multithread=TRUE, verbose = TRUE)

# Plot to inspect (should show black points over orange line)
plotErrors(errF.P2, nominalQ=TRUE)
plotErrors(errR.P2, nominalQ=TRUE)

### Sample Inference (using the PRE-LEARNED error models)
dadaFs.P1 <- dada(derepFs.P1, err=errF.P1, multithread=TRUE)  # Uses errF from learnErrors
dadaRs.P1 <- dada(derepRs.P1, err=errR.P1, multithread=TRUE)  # Uses errR from learnErrors

dadaFs.P2 <- dada(derepFs.P2, err=errF.P2, multithread=TRUE)  
dadaRs.P2 <- dada(derepRs.P2, err=errR.P2, multithread=TRUE) 


### Sample Inference
dadaFs.P1 <- dada(derepFs.P1, err=errF.P1, multithread=TRUE)
dadaRs.P1 <- dada(derepRs.P1, err=errR.P1, multithread=TRUE)

dadaFs.P2 <- dada(derepFs.P2, err=errF.P2, multithread=TRUE)
dadaRs.P2 <- dada(derepRs.P2, err=errR.P2, multithread=TRUE)


### Merge Paired Reads
mergers.P1 <- mergePairs(dadaFs.P1, derepFs.P1, dadaRs.P1, derepRs.P1, minOverlap = 20, verbose=TRUE)

mergers.P2 <- mergePairs(dadaFs.P2, derepFs.P2, dadaRs.P2, derepRs.P2, minOverlap = 20, verbose=TRUE)


### Construct sequence table 
seqtab.P1 <- makeSequenceTable(mergers.P1)

seqtab.P2 <- makeSequenceTable(mergers.P2)

### Merge the seqtab tables
seqtab <- mergeSequenceTables(seqtab.P1, seqtab.P2, tryRC = TRUE)


# # # FOR 12S!! ONLY - additional trimming step
obect <- as.data.frame(colnames(seqtab)) %>%
  dplyr::rename("sequence" = 1) %>%
  mutate(seq_length = nchar(sequence)) %>%
  mutate(new_sequence = substr(sequence, 1, 185)) %>%
  mutate(seq2_length = nchar(new_sequence))

colnames(seqtab) <- obect$new_sequence

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track.P1 <- cbind(out.P1, sapply(dadaFs.P1, getN), sapply(dadaRs.P1, getN), sapply(mergers.P1, getN), rowSums(seqtab.nochim[ sample.names.P1, , drop = FALSE ]))
colnames(track.P1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.P1) <- sample.names.P1

getN <- function(x) sum(getUniques(x))
track.P2 <- cbind(out.P2, sapply(dadaFs.P2, getN), sapply(dadaRs.P2, getN), sapply(mergers.P2, getN), rowSums(seqtab.nochim[ sample.names.P2, , drop = FALSE ]))
colnames(track.P2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.P2) <- sample.names.P2

### Appends the two track matrices
track <- rbind(track.P1, track.P2)

### Assign Taxonomy
taxam <- assignTaxonomy(seqtab.nochim, taxref, tryRC = TRUE, minBoot = 95)

# Assign Species
genus.species <- assignSpecies(seqtab.nochim, speciesref)

# Ensures all taxonomic levels agree for all rows within an species assignment
taxa <- as.data.frame(taxam)

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
save(seqtab.nochim, freq.nochim, track, out, taxa, file = "WADE003-arcticpred_dada2_QAQC_12S_output-batched.Rdata")

getwd()





#load("WADE003-arcticpred_dada2_QAQC_12S_output.Rdata")

### Follow up on unassigned sequences-------------------------------------------------------------------------
start_blast <- Sys.time()

# adds API key into environment
rentrez::set_entrez_key("53ea9c3064d51857214d415c23cee00a6508")

# Find sequences with NA below kingdom
unassigned_idx <- which(is.na(taxa[, "Species"]))  # or any level you care about

# Extract sequence names
seqs_unassigned <- colnames(seqtab.nochim)[unassigned_idx]

# Convert to DNAStringSet for export
seqs_unassigned_fasta <- DNAStringSet(seqs_unassigned)
names(seqs_unassigned_fasta) <- seqs_unassigned  # keep sequences as IDs

writeXStringSet(seqs_unassigned_fasta, "/gscratch/coenv/mball3/WADE003-arctic-pred/unassigned_seqs.fasta")

###uploaded this file to BLAST, then downloaded hits as a csv

unassigned_hits <- read_csv("/gscratch/coenv/mball3/WADE003-arctic-pred/unassigned_ALL_BLAST.csv", col_names = FALSE)

accessions <- unique(unassigned_hits$X2)

safe_esummary <- purrr::safely(function(acc) {
  entrez_summary(db = "nuccore", id = acc)
})

tax_summaries <- map(accessions, function(acc) {
  res <- safe_esummary(acc)
  if (!is.null(res$error)) {
    # optional: message about the bad accession
    message("Failed for accession: ", acc, " | error: ", res$error)
    return(NULL)
  }
  s <- res$result
  tibble(
    accession = acc,
    title     = s$title %||% NA_character_,
    taxid     = s$taxid %||% NA_character_,
    organism  = s$organism %||% NA_character_
  )
})

tax_df <- bind_rows(tax_summaries)

# Re-run with correct column reference
unassigned_taxa <- unassigned_hits %>%
  left_join(tax_df, by = c("X2" = "accession")) %>%
  filter(!is.na(organism)) %>%
  group_by(`X1`) %>%  # Backticks for column name
  mutate(
    organism_clean = str_extract(organism, "^[A-Z][a-z]+"),
    Genus = ifelse(n_distinct(organism_clean, na.rm = TRUE) > 1,
                   organism_clean[which.max(as.numeric(`X3`))],
                   organism_clean[1]),
    Species = "sp."
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(`X1`, Genus, Species, organism, everything())

head(unassigned_taxa)

# ####Look closely at this table to make sure the hits make sense! Maybe before the slice_head step
# 
# #now get full taxonomy

# get unique, non-NA taxids from your tax_df
taxids <- unique(tax_df$taxid[!is.na(tax_df$taxid)])
taxids <- as.character(taxids)

length(taxids)

# simple progress bar + fetch loop
pb <- progress_bar$new(total = length(taxids))
lineage_list <- lapply(taxids, function(tid) {
  pb$tick()
  Sys.sleep(0.3)
  tryCatch(
    entrez_fetch(db = "taxonomy", id = tid, rettype = "xml"),
    error = function(e) NA_character_
  )
})
pb$terminate()

valid_lineage <- lineage_list[!sapply(lineage_list, is.null) & !sapply(lineage_list, is.na)]

parse_lineage <- function(xml_string) {
  doc <- read_xml(xml_string)
  taxid <- xml_text(xml_find_first(doc, "//Taxon/TaxId"))
  lineage_nodes <- xml_find_all(doc, "//LineageEx/Taxon")
  ranks <- xml_text(xml_find_all(lineage_nodes, "Rank"))
  names <- xml_text(xml_find_all(lineage_nodes, "ScientificName"))
  lineage <- setNames(names, ranks)
  Species <- xml_text(xml_find_first(doc, "//Taxon/ScientificName"))
  lineage["species"] <- Species
  
  df <- as.data.frame(as.list(lineage), stringsAsFactors = FALSE)
  df$taxid <- taxid
  df
}

pb2 <- progress_bar$new(total = length(valid_lineage))
taxonomy_df_list <- lapply(valid_lineage, function(xml) {
  pb2$tick()
  parse_lineage(xml)
})
pb2$terminate()

taxonomy_df <- bind_rows(taxonomy_df_list)

# Now impose canonical order and names for ranks
taxonomy_df <- taxonomy_df %>%
  select(domain, phylum, class, order, family, genus, species, taxid) %>%
  mutate(taxid = as.character(taxid)) %>%
  distinct() %>%
  dplyr::rename(
    Kingdom = domain,
    Phylum = phylum,
    Class = class,
    Order = order,
    Family = family,
    Genus = genus,
    Species = species
  )

# Make sure taxid types match before joining
unassigned_taxa <- unassigned_taxa %>%
  mutate(taxid = as.character(taxid))

new_taxa <- unassigned_taxa %>%
  select(-Genus, -Species) %>%
  left_join(taxonomy_df, by = "taxid") %>%
  select(X1, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames("X1")

# Convert taxa objects to data frames
taxa_df     <- as.data.frame(taxa, stringsAsFactors = FALSE)
new_taxa_df <- as.data.frame(new_taxa, stringsAsFactors = FALSE)

# sanity checks
head(rownames(taxa_df))
head(rownames(new_taxa_df))
colnames(taxa_df)
colnames(new_taxa_df)

# figure out which rows in taxa have NAs in Species
na_species_idx <- which(is.na(taxa_df$Species))
na_sequences <- rownames(taxa_df)[na_species_idx]

length(na_sequences)

# use new_taxa to overwrite those NAs
seqs_with_new_taxa <- intersect(na_sequences, rownames(new_taxa_df))
length(seqs_with_new_taxa)

# replace in taxa
taxa_df[seqs_with_new_taxa, ] <- new_taxa_df[seqs_with_new_taxa, ]
compiled_taxa <- as.matrix(taxa_df)

# blast_finish <- Sys.time()

### Save data
save(seqtab.nochim, freq.nochim, track, taxa, compiled_taxa, new_taxa, file = "WADE003-arcticpred_dada2_QAQC_12S_output-batched.wBLASt.Rdata")

