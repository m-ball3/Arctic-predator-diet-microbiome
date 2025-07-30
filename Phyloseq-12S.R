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

# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-130trunc7.Rdata")

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-MFU_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata; filters out NA's (shipment 1)
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
ps.12s <- subset_taxa(ps.12s, Kingdom!="Bacteria")
#ps.12s <- prune_samples(sample_sums(ps.12s) > 0, ps.12s)
#ps.12s <- subset_taxa(ps.12s, !is.na(Species))

# Remove samples with total abundance == 0
ps.12s <- prune_samples(sample_sums(ps.12s) > 0, ps.12s)

# Saves phyloseq obj
saveRDS(ps.12s, "ps.12s")

# Plots stacked bar plot of abundance
plot_bar(ps.12s, fill="Species")

## MERGE TO SPECIES HERE (TAX GLOM)
ps.12s = tax_glom(ps.12s, "Species", NArm = FALSE)

# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.12s, fill="Species")

# Calculates relative abundance of each species 
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

# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot <- plot_bar(ps12s.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot

# Plots with ADFG IDs
sp.rel.plot +
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")

gen.rel.plot <- plot_bar(ps12s.rel, fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

gen.rel.plot + 
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")

fam.rel.plot <- plot_bar(ps12s.rel, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot 

fam.rel.plot + 
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")

# Facet wrapped by predator species
### I WANT BOXES AROUND THE DIFFERENT FACETS
faucet <- plot_bar(ps12s.rel, x="LabID", fill="Genus") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  scale_x_discrete(labels = label_map) +
  labs(x = "ADFG ID")+
  guides(fill = guide_legend(title = "Genus"))


faucet

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.12s))
colnames(otu.abs) <- as.data.frame(tax_table(ps.12s))$Species

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
write.csv(otu.abs, "ADFG_12s_absolute_speciesxsamples-trunc130.csv", row.names = FALSE)
write.csv(otu.prop, "ADFG_12s_relative_speciesxsamples-trunc130.csv", row.names = FALSE)





















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
