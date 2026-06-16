# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# THIS IS THE THIRD STEP AFTER DADA2
## rownames-match.r and replicates-contaminated.r must be run before this!
# ------------------------------------------------------------------


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

# ------------------------------------------------------------------
# Loads in Data
# ------------------------------------------------------------------
# Loads in output from DADA2 + filtered seqtab.nochim and samdf from rownames-match.r
load("./Scripts/ubiome/rdata/replicates-contaminated_ubiome.RData")

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.norep))
names(dna) <- taxa_names(ps.norep)
ps.norep <- merge_phyloseq(ps.norep, dna)
taxa_names(ps.norep) <- paste0("ASV", seq(ntaxa(ps.norep)))

nsamples(ps.norep)

# Filters out any Mammalia
ps.filt <- subset_taxa(ps.norep, Class!="Mammalia")

# Remove samples with total abundance == 0
ps.filt <- prune_samples(sample_sums(ps.filt) > 0, ps.filt)

## MERGE TO SPECIES HERE (TAX GLOM)
ps.filt = tax_glom(ps.filt, "Genus", NArm = FALSE)

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)
lowcount.minor <- genefilter_sample(ps.filt, f1, A=1)
ps.minor <- prune_taxa(lowcount.minor, ps.filt)

# Filtering to remove taxa with less than 1% of diet in 4 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.filt, f1, A=2) # errors occur at A>2
ps.major <- prune_taxa(lowcount.filt, ps.filt)

# Plots stacked bar plot of abundance - to confirm presence of NA's
# plot_bar(ps.16s, fill="Species")

# Plots stacked bar plot of abundance - to confirm presence of NA's
abs <- plot_bar(ps.filt, fill="Class")
abs

# Saves absolute abundance plot
ggsave("Deliverables/ubiome/ABSOLUTE.png", plot = abs, width = 30, height = 8, units = "in", dpi = 300)


# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
ps.rel <- transform_sample_counts(ps.filt, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps.minor.rel <- transform_sample_counts(ps.minor, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps.major.rel <- transform_sample_counts(ps.major, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps.minor.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps.major.rel))), arr.ind = TRUE)

#SAVES
### needs updating from SRKW to ADFG ubiome!
# save(ps.16s, ps16s.rel, samdf_filt, seqtab.nochim_filt, taxam, track, out, freq.nochim, file = "srkw-ps.16s.RData")
# save(seqtab.nochim_filt, freq.nochim, track, taxam, ps16s.rel, ps.16s, file = "SRKW-diet-16SALL.Rdata")
# save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.filt, ps.16s.filt.rel, file = "SRKW-diet-FILT-16SALL.Rdata")
# save(seqtab.nochim_filt, freq.nochim, track, taxam, ps.16s.major, ps.16s.major.rel, file = "SRKW-diet-MAJOR-16SALL.Rdata")

# ------------------------------------------------------------------
# PLOTS RELATIVE ABUNDANCE
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
rel.plot <- plot_bar(ps.major.rel, fill="Class")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
rel.plot

pred.facet <- plot_bar(ps.major.rel, x = "LabID", fill = "Class") +
  facet_wrap(~Predator, ncol = 4, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines")
  )
pred.facet


#saves plots 
ggsave("Deliverables/ubiome/ubiome-majorclass.png", plot = rel.plot, width = 30, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/ubiome/ubiome-majorpredfacet.png", plot = pred.facet, width = 30, height = 8, units = "in", dpi = 300)


# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.major))
colnames(otu.abs) <- as.data.frame(tax_table(ps.major))$Genus

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf$Specimen.ID # error here

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps.major.rel))
colnames(otu.prop) <- as.data.frame(tax_table(ps.major.rel))$Genus

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf_filt$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV
write.csv(otu.abs, "./Deliverables/ALL/SRKW_absolute_speciesxsamples-MAJOR.csv", row.names = TRUE)
write.csv(otu.prop, "./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv", row.names = TRUE)


# Calculates the total percent abundance for each prey species (across all samples)
options(scipen=999)#turns off scientific notation

# Remove the first column if it's sample IDs
df_no_id <- otu.abs[ , -1]

# Remove any rows/columns with all zeros or NAs
df_no_id <- df_no_id[rowSums(df_no_id) > 0, colSums(df_no_id) > 0]

# Calculate total abundance (sum across all samples per species)
total_abundance <- colSums(df_no_id)
total_abundance_df <- as.data.frame(total_abundance)
sum <- sum(total_abundance_df)

23348669/29075704

# Convert to percent abundance
percent_abundance <- total_abundance / sum(total_abundance_df$total_abundance)

# View result
percent_abundance_df <- as.data.frame(percent_abundance)

