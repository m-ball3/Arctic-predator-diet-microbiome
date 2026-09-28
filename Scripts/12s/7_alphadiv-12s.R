# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS ALPHA DIVERSITY METRICS ON 12s DATA
# THIS IS THE FIRST STATISTICAL INVESTIGATION 
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
# ------------------------------------------------------------------------------

# install.packages(
#   "rbiom",
#   type = "binary",
#   repos = "https://cran.rstudio.com"
# )

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(stringr)
library(phyloseq)
library(Hmisc)
library(rbiom)
library(mia) # via bioconductor

getwd()
load("./Scripts/ubiome/rdata/ubiome.ps-major.Rdata")

# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------
metadata <- data.frame(sample_data(ps.major))

class(metadata)

# ------------------------------------------------------------------------------
# Richness
# ------------------------------------------------------------------------------

# Robbins
## the ratio of singletons to total taxa


# Observed 
## the number of ASVs present

# gets the otu table
asv_mat <- as(otu_table(ps.major), "matrix")

# Number of observed ASVs per sample
n_asvs <- rowSums(asv_mat > 0)

# Optional: data frame with sample IDs
asv_per_sample <- data.frame(
  LabID = names(n_asvs),
  n_ASVs = unname(n_asvs),
  row.names = NULL
)

asv_per_sample
hist(asv_per_sample$n_ASVs)

# Adds the observed features to the metadata table by LabID
metadata <- as.data.frame(left_join(metadata, asv_per_sample, by= "LabID"))

# Box-and-whisker plot, with individual samples overlaid
obs_plot <- ggplot(metadata, aes(x = Predator, y = n_ASVs, fill = "Predator")) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(
    x = "Host",
    y = "Number of ASVs per sample",
    title = "observed ASV richness by host"
  ) +
  theme_classic() +
  theme(legend.position = "none")

obs_plot
