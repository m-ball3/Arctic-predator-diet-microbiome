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

devtools::install_github("benjjneb/dada2", force = TRUE, lib = .libPaths()[1])
BiocManager::install("dada2", lib = .libPaths()[1], force = TRUE)
BiocManager::install("S4Vectors")
install.packages("seqinr") 

library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)
library(dplyr)

diet.seqs.file <- "/gscratch/coenv/mball3/WADE003-arctic-pred/rawdata/ubiome"
taxref <- "/gscratch/coenv/mball3/530-meta-analysis/silva_nr99_v138.2_toGenus_trainset.fa.gz"
# speciesref <- "/gscratch/coenv/mball3/WADE003-arctic-pred/ADDSPECIES_MURI_MFU_07_2025.fasta"

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
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 40, truncLen=c(240, 200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose=TRUE)

### Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### Learn Error Rates

### proposed fix for binned error rates, but currently not working
# ### Learn Error Rates (BINNED QUALITY SCORES METHOD using pre-learned approach)
# novaBinnedErrfun <- makeBinnedQualErrfun(c(2, 12, 24, 38))
# 
# # Learn error rates FIRST from filtered reads (using binned error function)
# errF <- learnErrors(filtFs, errorEstimationFunction = novaBinnedErrfun, multithread=TRUE)
# errR <- learnErrors(filtRs, errorEstimationFunction = novaBinnedErrfun, multithread=TRUE)
# 
# # Plot to inspect (should show black points over orange line)
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)
# 
# ### Sample Inference (using the PRE-LEARNED error models)
# dadaFs <- dada(derepFs, err=errF, multithread=TRUE)  # Uses errF from learnErrors
# dadaRs <- dada(derepRs, err=errR, multithread=TRUE)  # Uses errR from learnErrors

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

# # # # FOR 12S!! ONLY
# obect <- as.data.frame(colnames(seqtab)) %>%
#   rename("sequence" = 1) %>%
#   mutate(seq_length = nchar(sequence)) %>%
#   mutate(new_sequence = substr(sequence, 1, 185)) %>%
#   mutate(seq2_length = nchar(new_sequence))

# colnames(seqtab) <- obect$new_sequence

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
# genus.species <- assignSpecies(seqtab.nochim, speciesref)

# # Ensures all taxonomic levels agree for all rows within an species assignment
# taxa <- as.data.frame(taxam) %>% # saves taxa as a dataframe
#   rownames_to_column("ASV") %>%  # makes ASV rownames to a column
#   filter(is.na(Species)) %>% # keeps only those in Species column that are NA
#   left_join(as.data.frame(genus.species) %>% # adds genus.species df to tax_table by ASV column
#               rownames_to_column("ASV"),  by = "ASV") %>%  
#   unite(col = addSpecies, Genus.y,Species.y, sep = " ") %>%  # combines genus.species df columns into one and renames column addspecies
#   ungroup() %>%  # ungroups by ASV
#   mutate(addSpecies = case_when(addSpecies == "NA NA"~NA, TRUE~addSpecies)) %>% # changes addSpecies column to R NA if NA NA characters and to what is in addSpecies column when there is a species
#   mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%  # changes NA to spp. in addSpecies column
#   mutate(.grp = ifelse(is.na(addSpecies),  paste0("NA_grp_", row_number()),  addSpecies)) %>%  # creates new column where NAs are specified by rowname (NA_grp_5)
#   group_by(.grp) %>%  # groups by column .grp
#   mutate(Class = if (length(unique(Class)) > 1) NA else Class) %>% # makes sure column and .grp agree; if not -> NA
#   mutate(Order = if (length(unique(Order)) > 1) NA else Order) %>% 
#   mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>% 
#   mutate(Genus.x = if (length(unique(Genus.x)) > 1) NA else Genus.x) %>% 
#   ungroup() %>%  # ungroups
#   dplyr::rename("Genus" = Genus.x, "Species" = addSpecies) %>%  # renames weird column names to names that make sense
#   select(-Species.x) %>%  # removes Species.x column
#   bind_rows(as.data.frame(taxa) %>% # makes taxa a df and adds its columns
#               rownames_to_column("ASV") %>%  # adds the ASVs as a column instead of rownames
#               filter(!is.na(Species))) %>%  # adds all rows back in 
#   mutate(.grp = ifelse(is.na(Species),  # changes .grp column to ???
#                        paste0("NA_grp_", row_number()),  Species)) %>%  
#   group_by(.grp) %>%  # groups by .grp
#   fill(Order, Family, Genus, .direction = "updown") %>% # fills in NAs with the previous or next non NA value (adds correct value where disagreements were)
#   ungroup() %>%  # ungroups
#   select(-.grp) %>%  # removes .grp column
#   column_to_rownames("ASV") %>%  # puts the columns ASV back to rownames
#   as.matrix() # transforms df back to matrix

### Save data
save(seqtab.nochim, freq.nochim, track, out, taxam, file = "WADE003-arcticpred_dada2_QAQC_ubiome_output.Rdata")

# ### Follow up on unassigned sequences
# start_blast <- Sys.time()
# 
# # Find sequences with NA below kingdom
# unassigned_idx <- which(is.na(taxa[, "Phylum"]))  # or any level you care about
# 
# # Extract sequence names
# seqs_unassigned <- colnames(seqtab.nochim)[unassigned_idx]
# 
# # Convert to DNAStringSet for export
# seqs_unassigned_fasta <- DNAStringSet(seqs_unassigned)
# names(seqs_unassigned_fasta) <- seqs_unassigned  # keep sequences as IDs
# 
# writeXStringSet(seqs_unassigned_fasta, "/gscratch/coenv/mball3/530-meta-analysis/unassigned_seqs.fasta")
# 
# ###uploaded this file to BLAST, then downloaded hits as a csv
# 
# unassigned_hits <- read_csv("/gscratch/coenv/mball3/530-meta-analysis/unassigned_ALL_BLAST.csv", col_names = FALSE)
# 
# accessions <- unique(unassigned_hits$X2)
# 
# safe_esummary <- purrr::safely(function(acc) {
#   entrez_summary(db = "nuccore", id = acc)
# })
# 
# tax_summaries <- map(accessions, function(acc) {
#   res <- safe_esummary(acc)
#   if (!is.null(res$error)) {
#     # optional: message about the bad accession
#     message("Failed for accession: ", acc, " | error: ", res$error)
#     return(NULL)
#   }
#   s <- res$result
#   tibble(
#     accession = acc,
#     title     = s$title %||% NA_character_,
#     taxid     = s$taxid %||% NA_character_,
#     organism  = s$organism %||% NA_character_
#   )
# })
# 
# tax_df <- bind_rows(tax_summaries)

# tax_summaries <- lapply(accessions, function(acc) {
#   summary <- entrez_summary(db="nuccore", id=acc)
#   list(
#     accession = acc,
#     title = summary$title,
#     taxid = summary$taxid,
#     organism = summary$organism
#   )
# })
# 
# tax_df <- do.call(rbind, lapply(tax_summaries, as.data.frame))

# unassigned_taxa <- unassigned_hits %>% 
#   left_join(tax_df, by = c("X2" = "accession")) %>% 
#   filter(!is.na(organism)) %>% 
#   group_by(X1) %>% 
#   mutate(
#     organism_clean = str_extract(organism, "^[A-Z][a-z]+"),  # First word = Genus
#     Genus = ifelse(n_distinct(organism_clean, na.rm = TRUE) > 1, 
#                    organism_clean[which.max(as.numeric(X3))], 
#                    organism_clean[1]),
#     Species = "sp."
#   ) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   select(X1, Genus, Species, organism, everything())
# 
# 
# 
# 
# ####Look closely at this table to make sure the hits make sense! Maybe before the slice_head step
# 
# #now get full taxonomy
# 
# 
# # Re-run with correct column reference
# unassigned_taxa <- unassigned_hits %>% 
#   left_join(tax_df, by = c("X2" = "accession")) %>% 
#   filter(!is.na(organism)) %>% 
#   group_by(`X1`) %>%  # Backticks for column name
#   mutate(
#     organism_clean = str_extract(organism, "^[A-Z][a-z]+"),
#     Genus = ifelse(n_distinct(organism_clean, na.rm = TRUE) > 1, 
#                    organism_clean[which.max(as.numeric(`X3`))], 
#                    organism_clean[1]),
#     Species = "sp."
#   ) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   select(`X1`, Genus, Species, organism, everything())
# 
# head(unassigned_taxa)
# 
# 
# # Chunk taxids for lineage lookup (avoid rate limit)
# taxids <- tax_df$taxid[!is.na(tax_df$taxid)]
# n_taxids <- length(taxids)
# batch_size <- 50
# n_batches <- ceiling(n_taxids / batch_size)
# 
# cat("Fetching lineage for", n_taxids, "taxids in", n_batches, "batches\n")
# 
# pb <- progress_bar$new(total = n_taxids)
# lineage_list <- list()
# 
# for(b in 1:n_batches) {
#   start_idx <- (b-1)*batch_size + 1
#   end_idx <- min(b*batch_size, n_taxids)
#   batch_taxids <- taxids[start_idx:end_idx]
#   
#   batch_lineage <- lapply(batch_taxids, function(tid) {
#     Sys.sleep(0.4)  # Rate limiting
#     pb$tick()
#     tryCatch({
#       entrez_fetch(db="taxonomy", id=tid, rettype="xml")
#     }, error = function(e) NA)
#   })
#   
#   lineage_list <- c(lineage_list, batch_lineage)
# }
# 
# pb$terminate()

### Save data
save(seqtab.nochim, freq.nochim, track, out, taxa, file = "WADE003-arcticpred_dada2_QAQC_ubiome_output.Rdata")

