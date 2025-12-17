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
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-regionalDB.Rdata")

# ------------------------------------------------------------------
# CREATES PHYLOSEQ
# ------------------------------------------------------------------

# Creates master phyloseq object
ps.12s <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
    sample_data(samdf), 
    tax_table(taxa)
  )

# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES
# ------------------------------------------------------------------

# Identifies rows in Specimen.ID that appear more than once (replicates)
duplicated_ids <- samdf$Specimen.ID[duplicated(samdf$Specimen.ID) | duplicated(samdf$Specimen.ID, fromLast = TRUE)]
unique_dup_ids <- unique(duplicated_ids) # "EB24PH075-S" = WADE 111 and WADE 124

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
ps.12s.replicates <- subset_samples(ps.12s, Specimen.ID %in% unique_dup_ids)

# Removes species assignments less than 100 reads for readability
ps.12s.replicates <- prune_taxa(taxa_sums(ps.12s.replicates) > 0, ps.12s.replicates)

# Create the stacked bar plot
plot_bar(ps.12s.replicates, x = "LabID", fill = "Species")

# # Ensures samples removed in filtering are removed from samdf
# replicate_to_remove <- "WADE-003-124"
# samdf <- samdf[!rownames(samdf) %in% replicate_to_remove, ]

# Recreates master phyloseq object
ps.12s <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(taxa)
)%>% 
  subset_samples(LabID != replicate_to_remove)
sample_names(ps.12s)

# ------------------------------------------------------------------
# CLEANS PHYLOSEQ
# ------------------------------------------------------------------

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.12s))
names(dna) <- taxa_names(ps.12s)
ps.raw <- merge_phyloseq(ps.12s, dna)
taxa_names(ps.12s) <- paste0("ASV", seq(ntaxa(ps.12s)))

nsamples(ps.12s)

# Saves phyloseq obj (RAW)
save(ps.12s_ci, ps.12s_bering, file = "MY_FILE_NAME.Rdata")

# Filters out anything not in Actinopteri
ps.12s <- subset_taxa(ps.12s, Class == "Actinopteri")
nsamples(ps.12s)

# Remove samples with total abundance < 100
ps.12s <- prune_samples(sample_sums(ps.12s) >= 100, ps.12s)
sample_sums(ps.12s)
nsamples(ps.12s)

## MERGE TO SPECIES HERE (TAX GLOM)
ps.12s = tax_glom(ps.12s, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)
lowcount.filt <- genefilter_sample(ps.12s, f1, A=1)
ps.12s <- prune_taxa(lowcount.filt, ps.12s) # WONDREING IF I SHOULD NOW RENAME PS.12S !! <- DO THIS

# Saves phyloseq obj (AFTER FILTERING)
save(ps.12s, "ps.12s.regions")

# Plots stacked bar plot of abundance
plot_bar(ps.12s, fill="Species")

# Plots stacked bar plot of abundance - to confirm presence of NA's
plot_bar(ps.12s, fill="Species")

# Transforms read counts to relative abundance of each species 
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

TEST <- as.data.frame(tax_table(ps.12s))
# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot <- plot_bar(ps12s.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot

# Plots with ADFG IDs
ADFG.sp<- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.sp

gen.rel.plot <- plot_bar(ps12s.rel, fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

ADFG.gen<- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.gen

fam.rel.plot <- plot_bar(ps12s.rel, fill="Family")+
  theme_minimal() +
  labs(y= "Proportion")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot 

ADFG.fam <- plot_bar(ps12s.rel, x = "Specimen.ID", fill="Family")+
  theme_minimal() +
  labs(y= "Proportion")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.fam

# Facet wrapped by predator species
faucet <- plot_bar(ps12s.rel, x="LabID", fill="Species") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 

ADFG.faucet <- plot_bar(ps12s.rel, x="Specimen.ID", fill="Species") +
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
ggsave("Deliverables/12S/regions/WADE labels/1S-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/ADFG-12S-species.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions/WADE labels/12S-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/ADFG-12S-genus.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions/WADE labels/12S-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/ADFG-12S-family.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions/WADE labels/12S-species-by-pred.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/ADFG-12S-species-by-pred.111125.png", plot = ADFG.faucet, width = 16, height = 8, units = "in", dpi = 300)


# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.12s))

# Changes NA.1 to it's corresponding ASV
taxa.names <- as.data.frame(tax_table(ps.12s)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species = case_when(is.na(Species)~ASV,
                               TRUE~Species)) %>% 
  pull(Species)

colnames(otu.abs) <- taxa.names

tax_table <- as.data.frame(tax_table(ps.12s))

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.abs$Specimen.ID <- samdf[rownames(otu.abs), "Specimen.ID"]

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
# Removes column for relative abundance calc
otu_counts <- otu.abs[, -1]

# row-wise proportions
otu.prop <- otu_counts / rowSums(otu_counts)

# add Specimen.ID back and reorder columns
otu.prop$Specimen.ID <- otu.abs$Specimen.ID
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# replace NaN (rows that were all zero) with 0
otu.prop[is.na(otu.prop)] <- 0

# round numeric columns
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV and adds Lab ID to otu.abs
write.csv(otu.abs %>% 
            rownames_to_column("LabID"), "./Deliverables/12S/regions/ADFG_12s_absolute_speciesxsamples.csv", row.names = FALSE)


write.csv(otu.prop%>% 
            rownames_to_column("LabID"), "./Deliverables/12S/regions/ADFG_12s_relative_speciesxsamples-trunc130-4.csv", row.names = FALSE)


















# OLD CODE

# # Creates a dataframe that maps ADFG IDs to sequences in taxa
# 
# ## Converts seqtab.nochim to long format (SampleID = WADE sample ID, Sequence = sequence column names)
# seqtab_long <- as.data.frame(seqtab.nochim)
# seqtab_long$WADE_ID <- rownames(seqtab_long)
# 
# seqtab_long <- seqtab_long %>%
#   pivot_longer(
#     cols = -WADE_ID,
#     names_to = "Sequence",
#     values_to = "Abundance"
#   ) %>%
#   filter(Abundance > 0)   # Keep only entries where the sequence is present in the sample
# 
# ## Keeps only the Sequence and WADE_ID columns
# mapped.sequences <- seqtab_long[, c("Sequence", "WADE_ID")]
# 
# # Turns rownames into a column for joining to mapped.sequences
# mapped.sequences <- mapped.sequences %>%
#   left_join(
#     samdf %>% rownames_to_column("WADE_ID") %>% select(WADE_ID, Specimen.ID),
#     by = "WADE_ID"
#   )
# 
# mapped.taxa <- as.data.frame(merged.taxa)
# 
# # Add Sequence as a column from the row names
# mapped.taxa$Sequence <- rownames(mapped.taxa)
# 
# # Merge with mapped.taxa to add ADFG_ID for each Sequence
# mapped.taxa <- merge(mapped.taxa, mapped.sequences[, c("Sequence", "Specimen.ID")], by = "Sequence", all.x = TRUE)
# 
# # Reorders columns to put Sequence and ADFG_ID first
# other_cols <- setdiff(names(mapped.taxa), c("Sequence", "Specimen.ID"))
# mapped.taxa <- mapped.taxa[, c("Sequence", "Specimen.ID", other_cols)]
