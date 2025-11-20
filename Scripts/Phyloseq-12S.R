# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Sets up the Environment and Loads in data
# ------------------------------------------------------------------

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(dplyr)
library(dada2)
library(patchwork)

# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output.Rdata")


# ------------------------------------------------------------------
# FORMATS METADATASHEET FOR PHYLOSEQ OBJ
# ------------------------------------------------------------------
# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-MFU_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata; filters out NA's (shipment 1)
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!is.na(LabID))

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator" & makes all lowercase
samdf <- dplyr::rename(samdf, Predator = Species)
samdf$Predator <- tolower(samdf$Predator)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# Remove rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ]

# Then set row names to LabID
rownames(samdf) <- samdf$LabID

# Only keeps rows that appear in both metadata and seq.tab 
## AKA only samples that made it through all steps 
common_ids <- intersect(rownames(samdf), rownames(seqtab.nochim))
samdf <- samdf[common_ids, ]
seqtab.nochim <- seqtab.nochim[common_ids, ]

# Checks for duplicates and identical sample rownames in both

### should return FALSE
any(duplicated(rownames(samdf)))
any(duplicated(rownames(seqtab.nochim)))

### should return TRUE
all(rownames(samdf) %in% rownames(seqtab.nochim))
all(rownames(seqtab.nochim) %in% rownames(samdf))

# Checks for samples in metadata but not in OTU table
setdiff(rownames(samdf), rownames(seqtab.nochim))

# Checks for samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf))

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)


# ------------------------------------------------------------------
# CREATES PHYLOSEQ
# ------------------------------------------------------------------

# Creates master phyloseq object
ps.12s <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
    sample_data(samdf), 
    tax_table(merged.taxa)
  )

# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES
# ------------------------------------------------------------------

# Identifies rows in Specimen.ID that appear more than once (replicates)
duplicated_ids <- samdf$Specimen.ID[duplicated(samdf$Specimen.ID) | duplicated(samdf$Specimen.ID, fromLast = TRUE)]
unique_dup_ids <- unique(duplicated_ids) # "EB24PH075-S" = WADE 111 and WADE 124

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
ps.12s.replicates <- subset_samples(ps.12s, Specimen.ID %in% unique_dup_ids)

# Removes species assignments less than 100 reads for readability
ps.12s.replicates <- prune_taxa(taxa_sums(ps.12s.replicates) > 100, ps.12s.replicates)

# Create the stacked bar plot
plot_bar(ps.12s.replicates, x = "LabID", fill = "Species.y")

# Ensures samples removed in filtering are removed from samdf
replicate_to_remove <- "WADE-003-124"
samdf <- samdf[!rownames(samdf) %in% replicate_to_remove, ]

# RECREATES PHYLOSEQ OBJECT WITHOUT REPLICATES
# Creates master phyloseq object
ps.12s <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(merged.taxa)
)

# ------------------------------------------------------------------
# CLEANS PHYLOSEQ
# ------------------------------------------------------------------

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.12s))
names(dna) <- taxa_names(ps.12s)
ps.raw <- merge_phyloseq(ps.12s, dna)
taxa_names(ps.12s) <- paste0("ASV", seq(ntaxa(ps.12s)))

nsamples(ps.12s)

# Filters out anything not in Actinopteri
ps.12s <- subset_taxa(ps.12s, Class == "Actinopteri")
nsamples(ps.12s)

# Remove samples with total abundance < 100
ps.12s <- prune_samples(sample_sums(ps.12s) >= 100, ps.12s)
sample_sums(ps.12s)
nsamples(ps.12s)

# Saves phyloseq obj
saveRDS(ps.12s, "ps.12s")

# Plots stacked bar plot of abundance
plot_bar(ps.12s, fill="Species.y")

## MERGE TO SPECIES HERE (TAX GLOM)
ps.12s = tax_glom(ps.12s, "Species.y", NArm = FALSE)

# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.12s, fill="Species.y")

# Transforms read counts to relative abundance of each species 
## Transforms NaN (0/0) to 0
ps12s.rel <- transform_sample_counts(ps.12s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps12s.rel))), arr.ind = TRUE)

# Creates a label map (WADE ID = ADFG ID)
label_map <- sample_data(ps12s.rel)$Specimen.ID
names(label_map) <- rownames(sample_data(ps12s.rel))

TEST <- as.data.frame(tax_table(ps.12s))
# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot <- plot_bar(ps12s.rel, fill="Species.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot

# Plots with ADFG IDs
ADFG.sp<- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Species.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.sp

gen.rel.plot <- plot_bar(ps12s.rel, fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

ADFG.gen<- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.gen

fam.rel.plot <- plot_bar(ps12s.rel, fill="Family")+
  theme_minimal() +
  labs(y= "Proportion")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot 

ADFG.fam <- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Family")+
  theme_minimal() +
  labs(y= "Proportion")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.fam

# Facet wrapped by predator species
faucet <- plot_bar(ps12s.rel, x="LabID", fill="Species.y") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 

ADFG.faucet <- plot_bar(ps12s.rel, x="Specimen.ID", fill="Species.y") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 

ADFG.faucet

#saves plots
ggsave("Deliverables/12S/WADE labels/1S-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/ADFG-12S-species.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/WADE labels/12S-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/ADFG-12S-genus.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/WADE labels/12S-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/ADFG-12S-family.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/WADE labels/12S-species-by-pred.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/ADFG-12S-species-by-pred.111125.png", plot = ADFG.faucet, width = 16, height = 8, units = "in", dpi = 300)


# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.12s))

# Changes NA.1 to it's corresponding ASV
taxa.names <- as.data.frame(tax_table(ps.12s)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species.y = case_when(is.na(Species.y)~ASV,
                               TRUE~Species.y)) %>% 
  pull(Species.y)

colnames(otu.abs) <- taxa.names

tax_table <- as.data.frame(tax_table(ps.12s))

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.abs$Specimen.ID <- samdf[rownames(otu.abs), "Specimen.ID"]

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps12s.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps12s.rel))$Species.y

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf[rownames(otu.prop), "Specimen.ID"]

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# # Replace NA column names in abs & prop
# na_cols <- which(is.na(colnames(otu.abs)))
# if(length(na_cols) > 0) colnames(otu.abs)[na_cols] <- "UNASSIGNED"

na_cols <- which(is.na(colnames(otu.prop)))
if(length(na_cols) > 0) colnames(otu.prop)[na_cols] <- "UNASSIGNED"

# Function to merge species variant columns (except for NA variants)
collapse_species <- function(df) {
  # Extract column names (assuming first is Specimen.ID)
  sp_cols <- colnames(df)[-1]
  # Only collapse if not 'NA' or an NA variant
  collapse_name <- function(x) {
    if (grepl('^NA(\\.|$)', x)) {
      return(x)
    } else {
      sub('\\.\\d+$', '', x)  # Remove trailing .number
    }
  }
  collapsed_names <- sapply(sp_cols, collapse_name)
  colnames(df)[-1] <- collapsed_names
  # Gather to long format
  df_long <- df %>% pivot_longer(-Specimen.ID, names_to = "Species", values_to = "Abundance")
  # Sum across variants for each specimen/species, keeping NA/NA.n separate
  df_sum <- df_long %>% group_by(Specimen.ID, Species) %>% summarize(Abundance = sum(Abundance), .groups = "drop")
  # Wide format back (Specimen.ID first)
  df_wide <- df_sum %>% pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)
  return(df_wide)
}

# Apply
otu.abs <- collapse_species(otu.abs)
otu.prop <- collapse_species(otu.prop)

# List of species names/keywords to remove
species_remove <- c(
  "delphinapterus", 
  "homo"
)

remove_species_cols <- function(df, remove_terms) {
  # TRUE if column matches ANY of the patterns in remove_terms
  match_any <- Reduce(`|`, lapply(remove_terms, function(term) grepl(term, colnames(df)[-1], ignore.case = TRUE)))
  # These columns to keep (not matched) plus 'Specimen.ID'
  keep_cols <- c('Specimen.ID', colnames(df)[-1][!match_any])
  df[, keep_cols, drop = FALSE]
}

# Apply for both 
otu.abs <- remove_species_cols(otu.abs, species_remove)
otu.prop <- remove_species_cols(otu.prop, species_remove)


# Writes to CSV
write.csv(otu.abs, "ADFG_12s_absolute_speciesxsamples-trunc130-4.csv", row.names = FALSE)
write.csv(otu.prop, "ADFG_12s_relative_speciesxsamples-trunc130-4.csv", row.names = FALSE)



















# OLD CODE

# # Creates a dataframe that maps ADFG IDs to sequences in taxa
# 
# ## Converts seqtab.nochim to long format (SampleID = WADE sample ID, Sequence = sequence column names)
# seqtab_long <- as.data.frame(seqtab.nochim)
# seqtab_long$WADE_ID <- rownames(seqtab_long)
# 
# seqtab_long <- seqtab_long %>%
#   pivot_longer(
#     cols = -WADE_ID,
#     names_to = "Sequence",
#     values_to = "Abundance"
#   ) %>%
#   filter(Abundance > 0)   # Keep only entries where the sequence is present in the sample
# 
# ## Keeps only the Sequence and WADE_ID columns
# mapped.sequences <- seqtab_long[, c("Sequence", "WADE_ID")]
# 
# # Turns rownames into a column for joining to mapped.sequences
# mapped.sequences <- mapped.sequences %>%
#   left_join(
#     samdf %>% rownames_to_column("WADE_ID") %>% select(WADE_ID, Specimen.ID),
#     by = "WADE_ID"
#   )
# 
# mapped.taxa <- as.data.frame(merged.taxa)
# 
# # Add Sequence as a column from the row names
# mapped.taxa$Sequence <- rownames(mapped.taxa)
# 
# # Merge with mapped.taxa to add ADFG_ID for each Sequence
# mapped.taxa <- merge(mapped.taxa, mapped.sequences[, c("Sequence", "Specimen.ID")], by = "Sequence", all.x = TRUE)
# 
# # Reorders columns to put Sequence and ADFG_ID first
# other_cols <- setdiff(names(mapped.taxa), c("Sequence", "Specimen.ID"))
# mapped.taxa <- mapped.taxa[, c("Sequence", "Specimen.ID", other_cols)]
