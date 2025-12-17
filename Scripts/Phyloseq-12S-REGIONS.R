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
library(patchwork)

# Loads dada2 output
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-COOKINLET.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-SBERING.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-ARCTIC.Rdata")

# ------------------------------------------------------------------
# CREATES PHYLOSEQ
# ------------------------------------------------------------------

# Creates REGIONAL phyloseq objects

ps.12s.cook <- phyloseq(
    otu_table(cook.seqtab, taxa_are_rows=FALSE), 
    sample_data(samdf), 
    tax_table(taxa.cook)
  ) # HAS DB NAMES 

ps.12s.sbering <- phyloseq(
  otu_table(sbering.seqtab, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(taxa.sbering)
)

ps.12s.arctic <- phyloseq(
  otu_table(arctic.seqtab, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(taxa.arctic)
)

# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES
# ------------------------------------------------------------------

# Identifies rows in Specimen.ID that appear more than once (replicates)
duplicated_ids <- samdf$Specimen.ID[duplicated(samdf$Specimen.ID) | duplicated(samdf$Specimen.ID, fromLast = TRUE)]
unique_dup_ids <- unique(duplicated_ids) # "EB24PH075-S" = WADE 111 and WADE 124

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
### replicates are only found in arctic db
ps.12s.replicates.arctic <- subset_samples(ps.12s.arctic, Specimen.ID %in% unique_dup_ids)

# Removes any taxa that have zero total reads across all samples
ps.12s.replicates.arctic <- prune_taxa(taxa_sums(ps.12s.replicates.arctic) > 0, ps.12s.replicates.arctic)

# Create the stacked bar plot
plot_bar(ps.12s.replicates.arctic, x = "LabID", fill = "Species")

# Specifies the replicate to remove
replicate_to_remove <- "WADE-003-124"

# Recreates REGIONAL phyloseq objects without unwanted replicate
ps.12s.arctic <- phyloseq(
  otu_table(arctic.seqtab, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(taxa.arctic)
)%>% 
  subset_samples(LabID != replicate_to_remove)
sample_names(ps.12s.arctic)

# ------------------------------------------------------------------
# CLEANS PHYLOSEQ
# ------------------------------------------------------------------

### shorten ASV seq names, store sequences as reference
dna.cook <- Biostrings::DNAStringSet(taxa_names(ps.12s.cook))
names(dna.cook) <- taxa_names(ps.12s.cook)
ps.cook.raw <- merge_phyloseq(ps.12s.cook, dna.cook) ## CHECK WITH AMY ABOUT THIS
taxa_names(ps.12s.cook) <- paste0("ASV", seq(ntaxa(ps.12s.cook)))

dna.sbering <- Biostrings::DNAStringSet(taxa_names(ps.12s.sbering))
names(dna.sbering) <- taxa_names(ps.12s.sbering)
ps.raw.sbering <- merge_phyloseq(ps.12s.sbering, dna.sbering) 
taxa_names(ps.12s.sbering) <- paste0("ASV", seq(ntaxa(ps.12s.sbering))) 

dna.arctic <- Biostrings::DNAStringSet(taxa_names(ps.12s.arctic))
names(dna.arctic) <- taxa_names(ps.12s.arctic)
ps.raw.arctic <- merge_phyloseq(ps.12s.arctic, dna.arctic)
taxa_names(ps.12s.arctic) <- paste0("ASV", seq(ntaxa(ps.12s.arctic)))

# compares number of samples
nsamples(ps.12s.cook)
nsamples(ps.12s.sbering)
nsamples(ps.12s.arctic)
#73 total samps

# Saves phyloseq obj per region (RAW)
save(ps.12s.cook, ps.12s.sbering, ps.12s.arctic, file = "ps.12s.regions.raw.Rdata")


# Filters out anything not in Actinopteri
ps.12s.cook <- subset_taxa(ps.12s.cook, Class == "Actinopteri")
nsamples(ps.12s.cook)

ps.12s.sbering <- subset_taxa(ps.12s.sbering, Class == "Actinopteri")
nsamples(ps.12s.sbering)

ps.12s.arctic <- subset_taxa(ps.12s.arctic, Class == "Actinopteri")
nsamples(ps.12s.arctic)

# Remove samples with total abundance < 100
ps.12s.cook <- prune_samples(sample_sums(ps.12s.cook) >= 100, ps.12s.cook)
sample_sums(ps.12s.cook)
nsamples(ps.12s.cook)

ps.12s.sbering <- prune_samples(sample_sums(ps.12s.sbering) >= 100, ps.12s.sbering)
sample_sums(ps.12s.sbering)
nsamples(ps.12s.sbering)

ps.12s.arctic <- prune_samples(sample_sums(ps.12s.arctic) >= 100, ps.12s.arctic)
sample_sums(ps.12s.arctic)
nsamples(ps.12s.arctic)

## MERGE TO SPECIES HERE (TAX GLOM)
ps.12s.cook = tax_glom(ps.12s.cook, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

ps.12s.sbering = tax_glom(ps.12s.sbering, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

ps.12s.arctic = tax_glom(ps.12s.arctic, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)

lowcount.filt.cook <- genefilter_sample(ps.12s.cook, f1, A=1)
ps.12s.cook.filt <- prune_taxa(lowcount.filt.cook, ps.12s.cook)

lowcount.filt.sbering <- genefilter_sample(ps.12s.sbering, f1, A=1)
ps.12s.sbering.filt <- prune_taxa(lowcount.filt.sbering, ps.12s.sbering)

lowcount.filt.arctic <- genefilter_sample(ps.12s.arctic, f1, A=1)
ps.12s.arctic.filt <- prune_taxa(lowcount.filt.arctic, ps.12s.arctic)


# Explores ASV assignments
asv.cook.rows <- as.data.frame(tax_table(ps.12s.cook.filt)) # LOST THE DB NAMES SOMEWHERE
asv.bering.rows <- as.data.frame(tax_table(ps.12s.sbering.filt))
asv.arctic.rows <- as.data.frame(tax_table(ps.12s.arctic.filt))

# makes rownames a column called rn
asv.cook.rows$rn <- rownames(asv.cook.rows)
asv.bering.rows$rn <- rownames(asv.bering.rows)
asv.bering.rows$rn <- rownames(asv.bering.rows)

# ADD CODE HERE TO MERGE DF BY ASV NAMES!
## NEED TO FIRST FIX DB NA ISSUE


# Saves phyloseq obj
save(ps.12s.cook.filt, ps.12s.sbering.filt, ps.12s.arctic.filt, file = "ps.12s.regions.filt.Rdata")

# Plots stacked bar plot of abundance
plot_bar(ps.12s.cook.filt, fill="Species")
plot_bar(ps.12s.sbering.filt, fill="Species")
plot_bar(ps.12s.arctic.filt, fill="Species")

# Transforms read counts to relative abundance of each species 
## Transforms NaN (0/0) to 0
ps12s.cook.rel <- transform_sample_counts(ps.12s.cook.filt, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps12s.sbering.rel <- transform_sample_counts(ps.12s.sbering.filt, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

ps12s.arctic.rel <- transform_sample_counts(ps.12s.arctic.filt, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps12s.cook.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps12s.sbering.rel))), arr.ind = TRUE)
which(is.nan(as.matrix(otu_table(ps12s.arctic.rel))), arr.ind = TRUE)

# Creates a label map (WADE ID = ADFG ID)
label_map <- sample_data(ps12s.cook.rel)$Specimen.ID
names(label_map) <- rownames(sample_data(ps12s.cook.rel))

label_map <- sample_data(ps12s.sbering.rel)$Specimen.ID
names(label_map) <- rownames(sample_data(ps12s.sbering.rel))

label_map <- sample_data(ps12s.arctic.rel)$Specimen.ID
names(label_map) <- rownames(sample_data(ps12s.arctic.rel))

TEST <- as.data.frame(tax_table(ps12s.cook.rel))
# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.cook.rel.plot <- plot_bar(ps12s.cook.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.cook.rel.plot

sp.sbering.rel.plot <- plot_bar(ps12s.sbering.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.sbering.rel.plot

sp.arctic.rel.plot <- plot_bar(ps12s.arctic.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.arctic.rel.plot

# ------------------------------------------------------------------
# WAITING TO UPDATE THIS REGIONALLY UNTIL WE FIGURE OUT HOW TO MERGE ALL BACK - MEB 12/17
# ------------------------------------------------------------------

# # Plots with ADFG IDs
# ADFG.sp<- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Species")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# ADFG.sp
# 
# gen.rel.plot <- plot_bar(ps12s.rel, fill="Genus")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# gen.rel.plot
# 
# ADFG.gen<- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Genus")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# ADFG.gen
# 
# fam.rel.plot <- plot_bar(ps12s.rel, fill="Family")+
#   theme_minimal() +
#   labs(y= "Proportion")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# fam.rel.plot 
# 
# ADFG.fam <- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Family")+
#   theme_minimal() +
#   labs(y= "Proportion")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# ADFG.fam
# 
# # Facet wrapped by predator species
# faucet <- plot_bar(ps12s.rel, x="LabID", fill="Species") +
#   facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     panel.spacing = unit(0.5, "lines"),
#     axis.title.x = element_text(margin = margin(t = 10))
#   ) 
# 
# ADFG.faucet <- plot_bar(ps12s.rel, x="Specimen.ID", fill="Species") +
#   facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     panel.spacing = unit(0.5, "lines"),
#     axis.title.x = element_text(margin = margin(t = 10))
#   ) 
# 
# ADFG.faucet
# 
# #saves plots
# ggsave("Deliverables/12S/regions/WADE labels/1S-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
# ggsave("Deliverables/12S/regions/ADFG-12S-species.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)
# 
# ggsave("Deliverables/12S/regions/WADE labels/12S-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
# ggsave("Deliverables/12S/regions/ADFG-12S-genus.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)
# 
# ggsave("Deliverables/12S/regions/WADE labels/12S-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
# ggsave("Deliverables/12S/regions/ADFG-12S-family.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)
# 
# ggsave("Deliverables/12S/regions/WADE labels/12S-species-by-pred.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
# ggsave("Deliverables/12S/regions/ADFG-12S-species-by-pred.111125.png", plot = ADFG.faucet, width = 16, height = 8, units = "in", dpi = 300)
# 

# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.cook.abs <- as.data.frame(otu_table(ps.12s.cook.filt))
otu.sbering.abs <- as.data.frame(otu_table(ps.12s.sbering.filt))
otu.arctic.abs <- as.data.frame(otu_table(ps.12s.arctic.filt))

# Changes NA.1 to it's corresponding ASV
cooktaxa.names <- as.data.frame(tax_table(ps.12s.cook.filt)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species = case_when(is.na(Species)~ASV,
                               TRUE~Species)) %>% 
  pull(Species)

colnames(otu.cook.abs) <- cooktaxa.names

cooktax_table <- as.data.frame(tax_table(ps.12s.cook.filt)) %>%
  select(-DB)

sberingtaxa.names <- as.data.frame(tax_table(ps.12s.sbering.filt)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species = case_when(is.na(Species)~ASV,
                             TRUE~Species)) %>% 
  pull(Species)

colnames(otu.sbering.abs) <- sberingtaxa.names

sberingtax_table <- as.data.frame(tax_table(ps.12s.sbering.filt)) %>%
  select(-DB)

arctictaxa.names <- as.data.frame(tax_table(ps.12s.arctic.filt)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species = case_when(is.na(Species)~ASV,
                             TRUE~Species)) %>% 
  pull(Species)

colnames(otu.arctic.abs) <- arctictaxa.names

arctictax_table <- as.data.frame(tax_table(ps.12s.arctic.filt)) %>%
  select(-DB)

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.cook.abs$Specimen.ID <- samdf[rownames(otu.cook.abs), "Specimen.ID"]
otu.sbering.abs$Specimen.ID <- samdf[rownames(otu.sbering.abs), "Specimen.ID"]
otu.arctic.abs$Specimen.ID <- samdf[rownames(otu.arctic.abs), "Specimen.ID"]

## Moves ADFG_SampleID to the first column
otu.cook.abs <- otu.cook.abs[, c(ncol(otu.cook.abs), 1:(ncol(otu.cook.abs)-1))]
otu.sbering.abs <- otu.sbering.abs[, c(ncol(otu.sbering.abs), 1:(ncol(otu.sbering.abs)-1))]
otu.arctic.abs <- otu.arctic.abs[, c(ncol(otu.arctic.abs), 1:(ncol(otu.arctic.abs)-1))]

# ------------------------------------------------------------------
# WAITING TO UPDATE THIS REGIONALLY UNTIL WE FIGURE OUT HOW TO MERGE ALL BACK - MEB 12/17
# ------------------------------------------------------------------

# # CREATES RELATIVE SAMPLES X SPECIES TABLE
# # Removes column for relative abundance calc
# otu_counts <- otu.abs[, -1]
# 
# # row-wise proportions
# otu.prop <- otu_counts / rowSums(otu_counts)
# 
# # add Specimen.ID back and reorder columns
# otu.prop$Specimen.ID <- otu.abs$Specimen.ID
# otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]
# 
# # replace NaN (rows that were all zero) with 0
# otu.prop[is.na(otu.prop)] <- 0
# 
# # round numeric columns
# is.num <- sapply(otu.prop, is.numeric)
# otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)
# 
# # Writes to CSV and adds Lab ID to otu.abs
# write.csv(otu.abs %>% 
#             rownames_to_column("LabID"), "./Deliverables/12S/regions/ADFG_12s_absolute_speciesxsamples.csv", row.names = FALSE)
# 
# write.csv(otu.prop%>% 
#             rownames_to_column("LabID"), "./Deliverables/12S/regions/ADFG_12s_relative_speciesxsamples-trunc130-4.csv", row.names = FALSE)
# 
# write.csv(tax_table%>% 
#             rownames_to_column("ASV"), "./Deliverables/12S/regions/ADFG_12s_tax_table.csv", row.names = FALSE)
