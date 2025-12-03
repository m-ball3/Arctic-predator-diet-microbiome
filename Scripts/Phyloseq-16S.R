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
library(tibble)


# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP1+2.Rdata")

# ------------------------------------------------------------------
# FORMATS METADATASHEET FOR PHYLOSEQ OBJ
# ------------------------------------------------------------------
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
      dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
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

# Removes (worst) replicate
sample_to_remove <- "WADE-003-118-C"

# Remove from metadata and OTU table early
samdf <- samdf[!rownames(samdf) %in% sample_to_remove, ]
seqtab.nochim <- seqtab.nochim[!rownames(seqtab.nochim) %in% sample_to_remove, ]

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)


# ------------------------------------------------------------------
# CREATES PHYLOSEQ
# ------------------------------------------------------------------

# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(merged.taxa))

# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES
# ------------------------------------------------------------------

# Identifies rows in Specimen.ID that appear more than once (replicates)
duplicated_ids <- samdf$Specimen.ID[duplicated(samdf$Specimen.ID) | duplicated(samdf$Specimen.ID, fromLast = TRUE)]
unique_dup_ids <- unique(duplicated_ids) # "PH22SH036-S" = WADE 115 and WADE 123

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
ps.16s.replicates <- subset_samples(ps.16s, Specimen.ID %in% unique_dup_ids)

# Removes species assignments less than 100 reads for readability
ps.16s.replicates <- prune_taxa(taxa_sums(ps.16s.replicates) > 0, ps.16s.replicates)

# Create the stacked bar plot
plot_bar(ps.16s.replicates, x = "LabID", fill = "Species.y")

# Ensures samples removed in filtering are removed from samdf
replicate_to_remove <- "WADE-003-115"
samdf <- samdf[!rownames(samdf) %in% replicate_to_remove, ]

# RECREATES PHYLOSEQ OBJECT WITHOUT REPLICATES
# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(merged.taxa))

# ------------------------------------------------------------------
# CLEANS PHYLOSEQ
# ------------------------------------------------------------------

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.16s))
names(dna) <- taxa_names(ps.16s)
ps.raw <- merge_phyloseq(ps.16s, dna)
taxa_names(ps.16s) <- paste0("ASV", seq(ntaxa(ps.16s)))

nsamples(ps.16s)

# Saves phyloseq obj
saveRDS(ps.16s, "ps.16s.raw")

# Filters out anything not in Actinopteri
ps.16s <- subset_taxa(ps.16s, Class == "Actinopteri")
nsamples(ps.16s)

# Remove samples with total abundance < 100
ps.16s <- prune_samples(sample_sums(ps.16s) >= 100, ps.16s)
nsamples(ps.16s)

# Ensures samples removed in filtering are removed from samdf
row_to_remove <- "WADE-003-146"
samdf <- samdf[!rownames(samdf) %in% row_to_remove, ]

# Saves phyloseq obj
saveRDS(ps.16s, "ps.16s")

## Merges same species
ps.16s = tax_glom(ps.16s, "Species.y", NArm = FALSE)

# Plots stacked bar plot of absolute abundance
plot_bar(ps.16s, x="Specimen.ID", fill="Species.y")

# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.16s, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)

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
ADFG.sp<- plot_bar(ps16s.rel, x = "Specimen.ID", fill="Species.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.sp

gen.rel.plot <- plot_bar(ps16s.rel, fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

ADFG.gen <- plot_bar(ps16s.rel, x = 'Specimen.ID', fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.gen

fam.rel.plot <- plot_bar(ps16s.rel, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot

ADFG.fam <- plot_bar(ps16s.rel, "Specimen.ID", fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.fam

# Facet wrapped by predator species
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

ADFG.faucet <- plot_bar(ps16s.rel, x ="Specimen.ID", fill = "Species.y") +
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
ggsave("Deliverables/16S/WADE labels/16S-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-species.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/WADE labels/16S-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-genus.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/WADE labels/16S-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-family.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/WADE labels/16S-species-by-pred.111125.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-species-by-pred.111125.png", plot = ADFG.faucet, width = 16, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.16s))

# Changes NA.1 to it's corresponding ASV
taxa.names <- as.data.frame(tax_table(ps.16s)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species.y = case_when(is.na(Species.y)~ASV,
                               TRUE~Species.y)) %>% 
  pull(Species.y)

colnames(otu.abs) <- taxa.names

tax_table <- as.data.frame(tax_table(ps.16s))

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
otu.prop[is.na(otu.prop)] <- 0

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

