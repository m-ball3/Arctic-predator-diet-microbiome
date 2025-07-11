# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# ------------------------------------------------------------------

## Sets Working Directory
# setwd("C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE lab/Arctic Predator/DADA2/DADA2 Outputs")
#setwd("DADA2/DADA2 Outputs")

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


# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output2.Rdata")

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-16S_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!is.na(LabID))

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator"
samdf <- dplyr::rename(samdf, Predator = Species)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# RemoveS rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ]

# Then set row names to LabID
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

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)

# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.16s))
names(dna) <- taxa_names(ps.16s)
ps.raw <- merge_phyloseq(ps.16s, dna)
taxa_names(ps.16s) <- paste0("ASV", seq(ntaxa(ps.16s)))

nsamples(ps.16s)

# Filters out any Mammalia and NA
ps.16s <- subset_taxa(ps.16s, Class!="Mammalia")
#ps.16s <- prune_samples(sample_sums(ps.16s) > 0, ps.16s)
#ps.16s <- subset_taxa(ps.16s, !is.na(Species))


## MERGE TO SPECIES HERE (TAX GLOM)
ps.16s = tax_glom(ps.16s, "Species")

# Plots stacked bar plot of abundance
plot_bar(ps.16s, fill="Species")

otu <- as.data.frame(otu_table(ps.16s))
# Calculates relative abundance of each species 
ps16s.rel <- transform_sample_counts(ps.16s, function(x) x/sum(x))

# Creates bar plot of relative abundance
rel.plot <- plot_bar(ps16s.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Extracts the sample data as a data frame
ADFG_sample_df <- as.data.frame(sample_data(ps.16s))

# Ensures the order of ADFG IDs matches the sample order in the plot
adfg_ids <- ADFG_sample_df$Specimen.ID[match(rel.plot$data$Sample, rownames(sample_df))]

# Overrides the x-axis labels with ADFG Sample IDs
rel.plot + scale_x_discrete(labels = adfg_ids)

# Facet wrapped by predator species
### I WANT BOXES AROUND THE DIFFERENT FACETS
faucet.plot <- plot_bar(ps16s.rel, fill = "Species") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  guides(fill = guide_legend(title = "Species"))

faucet.plot


# CREATES SAMPLES X SPECIES TABLE 

## Extracts the raw OTU count matrix (samples x ASVs)
otu_mat_raw <- as.data.frame(otu_table(ps.16s))
colnames(otu_mat_raw) <- as.data.frame(tax_table(ps.16s))$Species

##TRANSFORM FOR RELATIVE COUNTS


# tells that the species names are what we want colnames in df
#col.names(otu_mat_raw) <- tax_table(ps_object)$Species



## Extracts the taxonomy table (ASVs x taxonomy)
tax_mat <- as(tax_table(ps.16s), "matrix") ##LEAVE AS DF

## Gets species names for each ASV
species_names <- tax_mat[, "Species"]
species_names[is.na(species_names)] <- "Unknown"  # Replace NA species with 'Unknown'

## Transposes OTU matrix to ASVs x samples
otu_mat_t <- t(otu_mat_raw)

## Creates a data frame with species as a column
otu_df_t <- as.data.frame(otu_mat_t)
otu_df_t$Species <- species_names

## Aggregates by species, summing across ASVs for each sample
species_counts_t <- otu_df_t %>%
  group_by(Species) %>%
  summarise(across(everything(), sum))

## Transposes back to samples x species
species_counts <- as.data.frame(t(as.matrix(species_counts_t[,-1])))
colnames(species_counts) <- species_counts_t$Species

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
df_sample <- as.data.frame(sample_data(ps.16s))
species_counts$ADFG_SampleID <- df_sample$Specimen.ID

## Sanity check
head(species_counts)

### ADFG SAMPLE ID TO LEFT
species_counts$ADFG_SampleID <- sample_df$Specimen.ID

## Moves ADFG_SampleID to the first column
species_counts <- species_counts[, c(ncol(species_counts), 1:(ncol(species_counts)-1))]

# Writes to CSV
write.csv(species_counts, "ADFG_16s_speciesxsamples.csv", row.names = FALSE)


