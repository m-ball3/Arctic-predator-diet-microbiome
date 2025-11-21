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

# Checks relativbe abundance is calculated correctly 
equal.one <- combined_df %>%
  group_by(Sample, Marker) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  print(n = 22)

# # Ensures normalization after combination
# combined_df <- combined_df %>%
#   group_by(Sample, Marker) %>%
#   mutate(Abundance = Abundance / sum(Abundance)) %>%
#   ungroup()

# Extracts the sample data as a data frame
ADFG_sample_df <- as.data.frame(sample_data(ps.12s))


# Plots comparison
rel.plot <- ggplot(combined_df, aes(x = Sample, y = Abundance, fill = Species.y)) +
  geom_col(position = "stack") +
  facet_wrap(. ~ Marker, ncol=1, strip.position = "right") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(legend.position= "bottom")

rel.plot

# Ensures the order of ADFG IDs matches the sample order in the plot
adfg_ids <- ADFG_sample_df$Specimen.ID[match(rel.plot$data$Sample, rownames(ADFG_sample_df))]

# Overrides the x-axis labels with ADFG Sample IDs
p1 <- rel.plot + scale_x_discrete(labels = adfg_ids)
p1


ggsave("Deliverables/Comparisons/12S+16S.species-by-pred.png", plot = p1, width = 16, height = 10, units = "in", dpi = 300)


# 
# 
# rel.plot.sidebyside <- ggplot(combined_df, aes(x = Sample, y = Abundance, fill = Species.y)) +
#   geom_col(aes(group = Marker), position = position_dodge(width = 1)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
#   theme(legend.position = "bottom")
# 
# combined_df$Sample_Marker <- interaction(combined_df$Sample, combined_df$Marker)
# 
# ## NOT READY
# rel.plot.sidebyside <-
#   ggplot(combined_df, aes(x = Sample_Marker, y = Abundance, fill = Species.y)) +
#   geom_col() +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   theme(legend.position = "bottom")
# 
# 
# rel.plot.sidebyside
# 
# # Overrides the x-axis labels with ADFG Sample IDs
# ADFG.sidebysde <- rel.plot.sidebyside + scale_x_discrete(labels = adfg_ids)
# 
# 
# # ------------------------------------------------------------------
# # NOT READY
# # -----------------------------------------------------------------
# library(ggplot2)
# library(dplyr)
# library(forcats)

# # Get the sample order as a factor (desired order)
# # (This should have all the sample names, ordered as you want them to appear on x-axis)
# ordered_samples <- unique(combined_df$Sample)
# 
# combined_df$Sample <- factor(combined_df$Sample, levels = ordered_samples)
# combined_df$Marker <- factor(combined_df$Marker, levels = c("12S", "16S"))
# 
# # Now create a grouping for Sample+Marker and adjust ordering for side-by-side effect
# # You want: x = Sample, fill = Species.y, group = Marker, color = Marker
# # position_dodge separates 12S and 16S for each Sample
# 
# rel.plot.grouped <- ggplot(combined_df, aes(x = Sample, y = Abundance, fill = Species.y)) +
#   geom_col(aes(group = Marker, color = Marker), position = position_dodge(width = 0.8), width = 0.7) +
#   facet_wrap(~ Marker, nrow = 1) +   # Optional: if you want a separate facet for each marker
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "bottom"
#   )
# 
# # If you want no facet, just have pairs side-by-side for each sample:
# rel.plot.grouped <- ggplot(combined_df, aes(x = Sample, y = Abundance, fill = Species.y, group = Marker)) +
#   geom_col(position = position_dodge(width = 0.7), width = 0.6) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
# 
# # Use Specimen.ID as x labels (once per pair)
# ADFG_sample_df <- as.data.frame(sample_data(ps.12s))
# adfg_ids <- ADFG_sample_df$Specimen.ID[match(levels(combined_df$Sample), rownames(ADFG_sample_df))]
# 
# rel.plot.grouped + scale_x_discrete(labels = adfg_ids)
# 
# 
# 



















# 
# 
# # ------------------------------------------------------------------
# # Creates .csv of absolute and proportional samplexspecies
# # ------------------------------------------------------------------
# 
# # CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
# otu.abs <- cbind(combined_df$Specimen.ID, combined_df$Abundance)
# colnames(otu.abs) <- combined_df$Species
# 
# ## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
# otu.abs$Specimen.ID <- samdf$Specimen.ID
# 
# ## Moves ADFG_SampleID to the first column
# otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]
# 
# # CREATES RELATIVE SAMPLES X SPECIES TABLE
# otu.prop <- as.data.frame(otu_table(ps12s.rel))
# colnames(otu.prop) <- as.data.frame(tax_table(ps12s.rel))$Species
# 
# ## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
# otu.prop$Specimen.ID <- samdf$Specimen.ID
# 
# ## Moves ADFG_SampleID to the first column
# otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]
# 
# # Changes NaN to 0
# otu.prop[is.na(otu.prop)] <- 0
# 
# # Rounds to three decimal places
# is.num <- sapply(otu.prop, is.numeric)
# otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)
# 
# # Writes to CSV
# write.csv(otu.abs, "ADFG_12s_absolute_speciesxsamples.csv", row.names = FALSE)
# write.csv(otu.prop, "ADFG_12s_relative_speciesxsamples.csv", row.names = FALSE)
# 
