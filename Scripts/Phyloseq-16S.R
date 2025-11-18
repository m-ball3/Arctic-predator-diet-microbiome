# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# ------------------------------------------------------------------

## Sets Working Directory
setwd("Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs")

## Sets up the Environment and Libraries

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(dplyr)
library(tibble)


# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP1+2.Rdata")

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- sub("^((WADE-003-\\d+|WADE-003-\\d+-C|WADE-003-\\d+-UC))_.*", "\\1", rownames(seqtab.nochim))
rownames(seqtab.nochim) <- gsub("-16S_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator"
samdf <- dplyr::rename(samdf, Predator = Species)
samdf$Predator <- tolower(samdf$Predator)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% 
      select(Specimen.ID, Repeat.or.New.Specimen., LabID),
      by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# Removes rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ]

# Sets row names to LabID
rownames(samdf) <- samdf$LabID

# Only keeps rows that appear in both metadata and seq.tab 
## AKA only samples that made it through all steps 
common_ids <- intersect(rownames(samdf), rownames(seqtab.nochim))
samdf <- samdf[common_ids, ]
seqtab.nochim <- seqtab.nochim[common_ids, ]

# Checks for identical sample rownames in both
any(duplicated(rownames(samdf)))
any(duplicated(rownames(seqtab.nochim)))

all(rownames(samdf) %in% rownames(seqtab.nochim))
all(rownames(seqtab.nochim) %in% rownames(samdf))

# Samples in metadata but not in OTU table
setdiff(rownames(samdf), rownames(seqtab.nochim))

# Samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf))

sample_to_remove <- "WADE-003-118-C"

# Remove from metadata and OTU table early
samdf <- samdf[!rownames(samdf) %in% sample_to_remove, ]
seqtab.nochim <- seqtab.nochim[!rownames(seqtab.nochim) %in% sample_to_remove, ]

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)

# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(merged.taxa))

# Creates a dataframe that maps ADFG IDs to sequences in taxa

## Converts seqtab.nochim to long format (SampleID = WADE sample ID, Sequence = sequence column names)
seqtab_long <- as.data.frame(seqtab.nochim)
seqtab_long$WADE_ID <- rownames(seqtab_long)

seqtab_long <- seqtab_long %>%
  pivot_longer(
    cols = -WADE_ID,
    names_to = "Sequence",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0)   # Keep only entries where the sequence is present in the sample

## Keeps only the Sequence and WADE_ID columns
mapped.sequences <- seqtab_long[, c("Sequence", "WADE_ID")]

# Turns rownames into a column for joining to mapped.sequences
mapped.sequences <- mapped.sequences %>%
  left_join(
    samdf %>% rownames_to_column("WADE_ID") %>% select(WADE_ID, Specimen.ID),
    by = "WADE_ID"
  )

mapped.taxa <- as.data.frame(merged.taxa)

# Add Sequence as a column from the row names
mapped.taxa$Sequence <- rownames(mapped.taxa)

# Merge with mapped.taxa to add ADFG_ID for each Sequence
mapped.taxa <- merge(mapped.taxa, mapped.sequences[, c("Sequence", "Specimen.ID")], by = "Sequence", all.x = TRUE)

# Reorders columns to put Sequence and ADFG_ID first
other_cols <- setdiff(names(mapped.taxa), c("Sequence", "Specimen.ID"))
mapped.taxa <- mapped.taxa[, c("Sequence", "Specimen.ID", other_cols)]


### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.16s))
names(dna) <- taxa_names(ps.16s)
ps.raw <- merge_phyloseq(ps.16s, dna)
taxa_names(ps.16s) <- paste0("ASV", seq(ntaxa(ps.16s)))

nsamples(ps.16s)

# Filters out any Mammalia and NA
ps.16s <- subset_taxa(ps.16s, Class!="Mammalia")
ps.16s <- subset_taxa(ps.16s, Kingdom!="Bacteria")
#ps.16s <- prune_samples(sample_sums(ps.16s) > 0, ps.16s)
#ps.16s <- subset_taxa(ps.16s, !is.na(Species))


#sample 146 removed; need to remove from samdf
row_to_remove <- "WADE-003-146"
samdf <- samdf[!rownames(samdf) %in% row_to_remove, ]

# Remove samples with total abundance == 0
ps.16s <- prune_samples(sample_sums(ps.16s) > 0, ps.16s)

# Saves phyloseq obj
saveRDS(ps.16s, "ps.16s")


# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.16s, fill="Species.y")

## MERGE TO SPECIES HERE
ps.16s = tax_glom(ps.16s, "Species.y", NArm = FALSE)

# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.16s, fill="Species.y")

# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)

# Creates a label map (WADE ID = ADFG ID)
label_map <- sample_data(ps16s.rel)$Specimen.ID
names(label_map) <- rownames(sample_data(ps16s.rel))

# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot <- plot_bar(ps16s.rel, fill="Species.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot

# Plots with ADFG IDs
ADFG.sp<- sp.rel.plot +
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")
ADFG.sp

gen.rel.plot <- plot_bar(ps16s.rel, fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

ADFG.gen<- gen.rel.plot + 
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")
ADFG.gen

fam.rel.plot <- plot_bar(ps16s.rel, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot

ADFG.fam <- fam.rel.plot + 
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")
ADFG.fam

# Facet wrapped by predator species
### I WANT BOXES AROUND THE DIFFERENT FACETS
faucet <- plot_bar(ps16s.rel, fill = "Species.y") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 

ADFG.faucet <- faucet+
    scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")+
  guides(fill = guide_legend(title = "Species"))


ADFG.faucet

#saves plots 
ggsave("Deliverables/16S/16S-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-species.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/16S-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-genus.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/16S-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-family.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/16S-species-by-pred.111125.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-species-by-pred.111125.png", plot = ADFG.faucet, width = 16, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 


otu.abs <- as.data.frame(otu_table(ps.16s))
colnames(otu.abs) <- as.data.frame(tax_table(ps.16s))$Species.y

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf[rownames(otu.abs), "Specimen.ID"]

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps16s.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps16s.rel))$Species.y

## Adds ADFG Sample ID as a column
otu.prop$Specimen.ID <- samdf[rownames(otu.prop), "Specimen.ID"]

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)


# Replace NA column names in abs & prop
na_cols <- which(is.na(colnames(otu.abs)))
if(length(na_cols) > 0) colnames(otu.abs)[na_cols] <- "UNASSIGNED"

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

#---- 1. Apply to 16S ----#
otu.abs <- collapse_species(otu.abs)
otu.prop <- collapse_species(otu.prop)

# List of species names/keywords to remove
species_remove <- c(
  "Helicobacter",
  "Erignathus",
  "Pusa hispida"
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
write.csv(otu.abs, "ADFG_16s_absolute_speciesxsamples.csv", row.names = FALSE)
write.csv(otu.prop, "ADFG_16s_relative_speciesxsamples.csv", row.names = FALSE)

