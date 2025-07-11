# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# ------------------------------------------------------------------

## Sets Working Directory
# setwd("C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE lab/Arctic Predator/DADA2/DADA2 Outputs")
setwd("C:/Users/Intern/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs")

## Sets up the Environment and Libraries

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(dplyr)


# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("C:/Users/Intern/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP2_output3.Rdata")

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-MFU_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata; filters out NA's (shipment 1)
labdf <- read.csv("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!is.na(LabID))

samdf <- read.csv("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator"
samdf <- dplyr::rename(samdf, Predator = Species)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# RemoveS rows where LabID is NA (because shipment 1 was bad & thus not extracted)
## look into why I needed to do that again??
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

# Creates master phyloseq object
ps.12s <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
    sample_data(samdf), 
    tax_table(taxa)
  )

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.12s))
names(dna) <- taxa_names(ps.12s)
ps.raw <- merge_phyloseq(ps.12s, dna)
taxa_names(ps.12s) <- paste0("ASV", seq(ntaxa(ps.12s)))

nsamples(ps.12s)

# Filters out any Mammalia and NA
ps.12s <- subset_taxa(ps.12s, Class!="Mammalia")
ps.12s <- prune_samples(sample_sums(ps.12s) > 0, ps.12s)
#ps.12s <- subset_taxa(ps.12s, !is.na(Species))

# Plots stacked bar plot of abundance
plot_bar(ps.12s, fill="Species")

# Calculates relative abundance of each species 
ps12s.rel <- transform_sample_counts(ps.12s, function(otu) otu/sum(otu))

# Creates bar plot of relative abundance
rel.plot <- plot_bar(ps12s.rel, fill="Species")

rel.plot

# Extracts the sample data as a data frame
ADFG_sample_df <- as.data.frame(sample_data(ps.12s))

# Ensures the order of ADFG IDs matches the sample order in the plot
adfg_ids <- ADFG_sample_df$Specimen.ID[match(rel.plot$data$Sample, rownames(ADFG_sample_df))]

# Overrides the x-axis labels with ADFG Sample IDs
rel.plot + scale_x_discrete(labels = adfg_ids)



# cREATES SAMPLES X SPECIES TABLE 

## Extracts the raw OTU count matrix (samples x ASVs)
otu_mat_raw <- as(otu_table(ps.12s), "matrix")

## Extracts the taxonomy table (ASVs x taxonomy)
tax_mat <- as(tax_table(ps.12s), "matrix")

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
df_sample <- as.data.frame(sample_data(ps.12s))
species_counts$ADFG_SampleID <- df_sample$Specimen.ID

## Sanity check
head(species_counts)

### ADFG SAMPLE ID TO LEFT
species_counts$ADFG_SampleID <- sample_df$Specimen.ID

## Moves ADFG_SampleID to the first column
species_counts <- species_counts[, c(ncol(species_counts), 1:(ncol(species_counts)-1))]

# Writes to CSV
write.csv(species_counts, "ADFG_12s_speciesxsamples.csv", row.names = FALSE)























# Creates stacked bar plot 
## proportion of each species 
### samples on x
ps12s.prop <- transform_sample_counts(ps.12s, function(OTU) OTU/sum(OTU))
ps12s.bar <- transform_sample_counts(ps.12s, function(OTU) OTU/sum(OTU))

# Checks what taxa are present after transformation
taxa_sums(ps12s.bar)

otu_mat <- as(otu_table(ps12s.bar), "matrix")
rowSums(is.na(otu_mat))  # Should be zero for all ASVs
rowSums(otu_mat)         # Should be >0 for ASVs with data


plot_bar(ps12s.prop, x="LabID", fill="Species") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


# Facet wrapped by predator species
### I WANT BOXES AROUND THE DIFFERENT FACETS
plot_bar(ps12s.bar, x="LabID", fill="Species") +
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


# Compares stomach to fecal
sample_data(ps12s.bar)$Stomach.Goo <- factor(
  sample_data(ps12s.bar)$Stomach.Goo,
  levels = c("No", "Yes"),
  labels = c("Fecal", "Stomach")
)


plot_bar(ps12s.bar.no.beluga, x = "LabID", fill = "Species") +
  facet_wrap(~ Stomach.Goo, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  guides(fill = guide_legend(title = "Species"))



# Converts to proportional abundance

sample_variables(ps.12s)
otu_mat <- as(otu_table(ps.12s), "matrix")
sample_totals <- sample_sums(ps.12s)

# Extract the species assignments for each ASV
species_assignments <- tax_table(ps.12s)[, "Species"]

# Convert to character vector for assignment
species_names <- as.character(species_assignments[colnames(otu_mat)])

# Replace column names in otu_mat
colnames(otu_mat) <- species_names
colnames(otu_mat) <- make.unique(species_names)

first_asv <- colnames(otu_mat)[1]
tax_table(ps.12s)[first_asv, "Species"]



ps12s.propa <- transform_sample_counts(ps12s.prop, function(x) x / sum(x) )

plot_bar(ps12s.propa, x="LabID", fill="Species") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

### export csv for ampbias correction
# readcount.table <- as.data.frame(otu_table(ps.raw))
# taxon.table <- as.data.frame(tax_table(ps.raw)) %>% rownames_to_column(var = "ASV")
# metadata.table <- samdf %>% rownames_to_column(var = "Sample")
# reference.table <- as.data.frame(refseq(ps.raw)) %>% 
#   rownames_to_column("ASVname")
