# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# ------------------------------------------------------------------

## Sets Working Directory
# setwd("C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE lab/Arctic Predator/DADA2/DADA2 Outputs")
setwd("Arctic-predator-diet-microbiome/")

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


# Loads in phyloseq objects
ps.12s <- readRDS("ps.12s")
ps.16s <- readRDS("ps.16s")

# Gets the samples
sam.12s <- sample_names(ps.12s)
sam.16s <- sample_names(ps.16s)

# Selects only the samples that are present in both
sam.both <- intersect(sam.16s, sam.12s)

# Keeps only the samples that are in both
ps.12s <- prune_samples(sam.both, ps.12s)
ps.16s <- prune_samples(sam.both, ps.16s)

# Gets the relative abundance pf each species
rel.12s <- transform_sample_counts(ps.12s, function(x) x / sum(x))
rel.16s <- transform_sample_counts(ps.16s, function(x) x / sum(x))

df_12s <- psmelt(rel.12s)
df_16s <- psmelt(rel.16s)
df_12s$Marker <- "12S"
df_16s$Marker <- "16S"

# Creates the combined df
combined_df <- bind_rows(df_12s, df_16s)

# Plots absolute comparison
rel.plot <- ggplot(combined_df, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_col(position = "stack") +
  facet_wrap(. ~ Marker, ncol=1, strip.position = "right") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

rel.plot
# Extracts the sample data as a data frame
ADFG_sample_df <- as.data.frame(sample_data(ps.12s))

# Ensures the order of ADFG IDs matches the sample order in the plot
adfg_ids <- ADFG_sample_df$Specimen.ID[match(rel.plot$data$Sample, rownames(ADFG_sample_df))]

# Overrides the x-axis labels with ADFG Sample IDs
rel.plot + scale_x_discrete(labels = adfg_ids)

# ------------------------------------------------------------------
# Creates .csv of absolute and proportional samplexspecies
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- cbind(combined_df$Specimen.ID, combined_df$Abundance)
colnames(otu.abs) <- combined_df$Species

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.abs$Specimen.ID <- samdf$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps12s.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps12s.rel))$Species

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV
write.csv(otu.abs, "ADFG_12s_absolute_speciesxsamples-trunc110.csv", row.names = FALSE)
write.csv(otu.prop, "ADFG_12s_relative_speciesxsamples-trunc110.csv", row.names = FALSE)

