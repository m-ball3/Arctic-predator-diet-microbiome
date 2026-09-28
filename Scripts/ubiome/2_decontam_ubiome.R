# ------------------------------------------------------------------------------------
# DECONTAMINATION
# THIS IS THE SECOND STEP AFTER DADA2
## rownames-match.r and replicates-contaminated.r must be run before this!
# ------------------------------------------------------------------------------------

# Sets up environment
library(tidyverse)
library(phyloseq)
library(decontam)
library(vegan)
library(Biostrings)
library(openxlsx)

# Loads in data
load("./Scripts/ubiome/rdata/rownames-match_ubiome.RData")
labdf <- read.csv("./metadata/ADFG_dDNA_labwork_metadata.csv")

# ------------------------------------------------------------------------------------
# Creates a series of tables to compare read counts and ASVs in negs vs all samples
# ------------------------------------------------------------------------------------

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.raw))
names(dna) <- taxa_names(ps.raw)
ps.raw <- merge_phyloseq(ps.raw, dna)
taxa_names(ps.raw) <- paste0("ASV", seq(phyloseq::ntaxa(ps.raw)))

# saves the full ASV OTU table
otu.abs <- as.data.frame(otu_table(ps.raw))

# saves just the negs in an AVS OTU table
otu.abs.negs <- otu.abs %>%
  filter(
    grepl("neg",  rownames(.), ignore.case = TRUE))

# saves a table that just has the ASVs across samples that are in the negs
otu.abs.asvs.in.negs <- otu.abs[, colSums(otu.abs.negs) > 0]

# gets the sequences
seqs <- refseq(ps.raw)

# makes sure the sequences and the ASV#'s will correspond correctly
identical(taxa_names(ps.raw), names(seqs))

# creates a table mapping the sequence to the ASV#
asv_sequence_table <- data.frame(
  ASV = taxa_names(ps.raw),
  Sequence = as.character(seqs),
  row.names = NULL
)

# adds taxonomy to the asv_sequence table
taxams <- taxam %>%
  as.data.frame()%>%
  tibble::rownames_to_column("Sequence")

asv_sequence_table <- asv_sequence_table %>%
  dplyr::left_join(taxams, by = "Sequence")

# changes the rowname to a column named Sample
otu.abs <- rownames_to_column(otu.abs, var = "Sample")
otu.abs.negs <- rownames_to_column(otu.abs.negs, var = "Sample")
otu.abs.asvs.in.negs <- rownames_to_column(otu.abs.asvs.in.negs, var = "Sample")

otu.abs.negs <- otu.abs.negs[, colSums(otu.abs.negs > 0) > 0] # filters out all 0 ASVs across negs

# saves as a multi-column xlsx
write.xlsx(
  list(
    otu.abs = otu.abs,
    otu.abs.negs = otu.abs.negs,
    otu.abs.asvs.in.negs = otu.abs.asvs.in.negs,
    asv_sequence_table = asv_sequence_table
  ), rownames = TRUE,
  file = "./Deliverables/ubiome/contam_otu_absolute_data.xlsx"
)

# ------------------------------------------------------------------------------
# Pulls out just the contaminated samples
# ------------------------------------------------------------------------------

# Subsets to just shipment 5 (contaminated extraction)
# samdf for negs and shipment 5
samdf_contaminated <- samdf %>%
  filter(
    grepl("neg",  rownames(.), ignore.case = TRUE) |
      grepl("5", Shipment, ignore.case = TRUE)
  )

# labdf for negs and shipment 5
labdf_contaminated <- labdf %>%
  filter(
  grepl("neg",  LabID, ignore.case = TRUE) |
    grepl("5", Shipment, ignore.case = TRUE)
)

# makes lab id the rownames
rownames(labdf_contaminated) <- labdf_contaminated$LabID

# makes the rows in samdf_contaminated match those in labdf_contaminated
contaminated_samps <- rownames(samdf_contaminated)

# Subset labdf_contaminated to matching rows
labdf_contaminated <- labdf_contaminated[
  rownames(labdf_contaminated) %in% contaminated_samps,
  ,
  drop = FALSE
]

# subset a phyloseq object to those samples

ps.contaminated <- prune_samples(contaminated_samps, ps.raw)

# ### shorten ASV seq names, store sequences as reference
# dna <- Biostrings::DNAStringSet(taxa_names(ps.contaminated))
# names(dna) <- taxa_names(ps.contaminated)
# ps.contaminated <- merge_phyloseq(ps.contaminated, dna)
# taxa_names(ps.contaminated) <- paste0("ASV", seq(ntaxa(ps.contaminated)))

phyloseq::nsamples(ps.contaminated)

## MERGE TO SPECIES HERE (TAX GLOM)
ps.contaminated.glom = tax_glom(ps.contaminated, "Genus", NArm = FALSE)

## Transforms NaN (0/0) to 0
ps.contaminated.prop <- transform_sample_counts(ps.contaminated.glom, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# PLOTS ------------------------------------------------------------------------
# Plots with WADE IDs
ps.contaminated.plot <- plot_bar(ps.contaminated, fill="Genus")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = 'right')
ps.contaminated.plot

# Plots with WADE IDs
rel.plot <- plot_bar(ps.contaminated.prop, fill="Class")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = 'right')
rel.plot

# ggsave("Deliverables/ubiome/ubiome-negmock-Class.png", plot = rel.plot, width = 30, height = 30, units = "in", dpi = 300)

# TABLES -----------------------------------------------------------------------
# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.contaminated))
colnames(otu.abs) <- as.data.frame(tax_table(ps.contaminated))$Family

# Checks for differences = WADE-003-178
setdiff(rownames(samdf_contaminated), rownames(otu.abs))
setdiff(rownames(otu.abs), rownames(samdf_contaminated))

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf_contaminated$Specimen.ID # error here

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps.contaminated.prop))
colnames(otu.prop) <- as.data.frame(tax_table(ps.contaminated.prop))$Family

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf_contaminated$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV
# write.csv(otu.abs, "./Deliverables/ubiome/ubiome_absolute_speciesxsamples-negmock-FAM.csv", row.names = TRUE)
# write.csv(otu.prop, "./Deliverables/ubiome/ubiome_relative_speciesxsamples-negmock-FAM.csv", row.names = TRUE)


# Creates a table that only keeps the ASVs that are in the negs 

neg.ASVs <- read.csv("./Deliverables/ubiome/ubiome_absolute_speciesxsamples-NEGS.csv")
all.ASVs <- read.csv("./Deliverables/ubiome/ubiome_absolute_speciesxsamples-NEGS.csv")

# Visualizes shipment 5 compared to the contaminated negatives

# ------------------------------------------------------------------------------------
# USES FOUR METHODS FOR DECONTAMINATION
# 1. manually deleting the contaminated reads
# 2. using R package decontam - prevalence method
# 3. using R package decontam - frequency method
# 4. using R package decontam - combination of prevalence & frequency methods
# ------------------------------------------------------------------------------------

# MANUALLY DELETES READS--------------------------------------------------------

# Extract the OTU table as a matrix (taxa x samples)
otu_mat <- as.matrix(otu_table(ps.contaminated))

# Ensure taxa are rows
if (!taxa_are_rows(ps.contaminated)) {
  otu_mat <- t(otu_mat)
}

#gets the negative samples 
neg_ids <- c("neg1", "neg2", "neg2-2")  

neg_ids <- intersect(neg_ids, colnames(otu_mat))

##  gets mean read count for each species in negs
neg_mat <- otu_mat[ , neg_ids, drop = FALSE]

# Mean read count per species (row) across the negatives
neg_means <- rowMeans(neg_mat)


##  uses that to delete the mean neg read count from each corresponding species in samps
# Subtract mean contamination from every sample
otu_decon <- sweep(otu_mat, 1, neg_means, FUN = "-")

# Floor negative values at 0
otu_decon[otu_decon < 0] <- 0

### WADE-003-151/EB10PH008 should now = 1265 for Clostridaceae

# Rebuild an otu_table with the same orientation as original
otu_decon_tab <- otu_table(otu_decon, taxa_are_rows = TRUE)
otu_decon_tab <- t(otu_decon_tab)

ps.decon <- phyloseq(
  otu_decon_tab,
  sample_data(ps.contaminated),
  tax_table(ps.contaminated)
)

# checks table
# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs.decon <- as.data.frame(otu_table(ps.decon))
colnames(otu.abs.decon) <- as.data.frame(tax_table(ps.decon))$Family

# Checks for differences = WADE-003-178
setdiff(rownames(samdf_contaminated), rownames(otu.abs.decon))
setdiff(rownames(otu.abs.decon), rownames(samdf_contaminated))

## Adds ADFG Sample ID as a column
otu.abs.decon$Specimen.ID <- samdf_contaminated$Specimen.ID # error here

## Moves ADFG_SampleID to the first column
otu.abs.decon <- otu.abs.decon[, c(ncol(otu.abs.decon), 1:(ncol(otu.abs.decon)-1))]

# DEALS WITH CONTAMINATION USING R PACKAGE DECONTAM PREVALENCE METHOD -----------------------------

# creates a copy of the contaminated samples ps to run decontam on
ps_for_decontam <- ps.contaminated

# Label negatives
sample_data(ps_for_decontam)$is_neg <- sample_names(ps_for_decontam) %in% neg_ids
table(sample_data(ps_for_decontam)$is_neg)  # 55 FALSE, 3 TRUE

# Run decontam on the phyloseq object (prevalence method)
contam_prev <- isContaminant(
  ps_for_decontam,
  method = "prevalence",
  neg    = "is_neg",
  threshold = 0.1
)

# sanity checks
table(contam_prev$contaminant)   # 2344 FALSE, 53 TRUE
length(contam_prev$contaminant)  # 2397
head(rownames(contam_prev))      # "ASV1", "ASV2", ...

# gets the contaminated ASVs
contam_asvs <- rownames(contam_prev)[contam_prev$contaminant]

# prunes contaminated ASVs from ps object and creates new ps
ps.decontam <- prune_taxa(
  !taxa_names(ps_for_decontam) %in% contam_asvs,
  ps_for_decontam
)

# sanity checks
phyloseq::ntaxa(ps_for_decontam)  # 2397
phyloseq::ntaxa(ps.decontam)      # should be 2397 - 53 = 2344

# checks table
# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs.decontam <- as.data.frame(otu_table(ps.decontam))
colnames(otu.abs.decontam) <- as.data.frame(tax_table(ps.decontam))$Family

# Checks for differences = WADE-003-178
setdiff(rownames(samdf_contaminated), rownames(otu.abs.decontam))
setdiff(rownames(otu.abs.decontam), rownames(samdf_contaminated))

## Adds ADFG Sample ID as a column
otu.abs.decontam$Specimen.ID <- samdf_contaminated$Specimen.ID # error here

## Moves ADFG_SampleID to the first column
otu.abs.decontam <- otu.abs.decontam[, c(ncol(otu.abs.decontam), 1:(ncol(otu.abs.decontam)-1))]


# DEALS WITH CONTAMINATION USING R PACKAGE DECONTAM FREQUENCY METHOD COMBINED WITH PREVALENCE METHOD -----------------------------

# Add DNA conc as a numeric column in sample_data
sample_data(ps_for_decontam)$Final_quant_ubiome <- labdf_contaminated$Final_quant_ubiome[
  match(sample_names(ps_for_decontam), rownames(labdf_contaminated))
]

# Run decontam on the phyloseq object (frequency method)
contam_freq <- isContaminant(
  ps_for_decontam,
  method = "frequency",
  conc   = "Final_quant_ubiome",
  threshold = 0.1
)

table(contam_freq$contaminant)

# Combined contaminant flag
combined_contam <- contam_prev$contaminant | contam_freq$contaminant

contam_asvs_combined <- rownames(contam_prev)[combined_contam]

ps.decontam_freandprev <- prune_taxa(
  !taxa_names(ps_for_decontam) %in% contam_asvs_combined,
  ps_for_decontam
)

test <- as.data.frame(combined_contam)


# DEALS WITH CONTAMINATION USING R PACKAGE 3 -----------------------------------

# STATS TO COMPARE -------------------------------------------------------------

# Extract matrices with samples as rows
otu.manual   <- as(otu_table(ps.decon), "matrix")
otu.decontam <- as(otu_table(ps.decontam), "matrix")
otu.contam   <- as(otu_table(ps.contaminated), "matrix")
otu.combined <- as(otu_table(ps.decontam_freandprev), "matrix")

# Align sample IDs
common_samps <- Reduce(intersect, list(
  rownames(otu.manual),
  rownames(otu.decontam),
  rownames(otu.contam), 
  rownames(otu.combined)
))

otu.manual   <- otu.manual[common_samps, , drop = FALSE]
otu.decontam <- otu.decontam[common_samps, , drop = FALSE]
otu.contam   <- otu.contam[common_samps, , drop = FALSE]
otu.combined <- otu.combined[common_samps, , drop = FALSE]

# Union of all ASV columns
common_taxa <- union(colnames(otu.manual), colnames(otu.decontam))
common_taxa <- union(common_taxa, colnames(otu.contam))
common_taxa <- union(common_taxa, colnames(otu.combined))

# Helper to expand a matrix to the full set of taxa
expand_to_common_taxa <- function(mat, common_taxa) {
  # Create an all-zero matrix with common_taxa columns
  expanded <- matrix(
    0,
    nrow = nrow(mat),
    ncol = length(common_taxa),
    dimnames = list(rownames(mat), common_taxa)
  )
  # Find taxa present in this matrix
  present_taxa <- intersect(common_taxa, colnames(mat))
  # Copy over the counts for those taxa
  expanded[ , present_taxa] <- mat[ , present_taxa, drop = FALSE]
  expanded
}

otu.manual   <- expand_to_common_taxa(otu.manual, common_taxa)
otu.decontam <- expand_to_common_taxa(otu.decontam, common_taxa)
otu.contam   <- expand_to_common_taxa(otu.contam, common_taxa)
otu.combined <- expand_to_common_taxa(otu.combined, common_taxa)


method.manual <- as.data.frame(otu.manual)
method.manual$Sample <- rownames(method.manual)
method.manual$method <- "manual"

method.decontam <- as.data.frame(otu.decontam)
method.decontam$Sample <- rownames(method.decontam)
method.decontam$method <- "decontam"

method.contam <- as.data.frame(otu.contam)
method.contam$Sample <- rownames(method.contam)
method.contam$method <- "contam"

method.combined <- as.data.frame(otu.combined)
method.combined$Sample <- rownames(method.combined)
method.combined$method <- "combined"


stacked <- rbind(method.manual, method.decontam, method.contam, method.combined)
abund <- as.matrix(stacked[ , !(names(stacked) %in% c("Sample", "method"))])

## ALPHA DIVERSITY -----------------------------------------------------------------------
shannon <- diversity(abund, index = "shannon", MARGIN = 1)  # base e by default

stacked$Shannon <- shannon

ggplot(stacked, aes(x = method, y = Shannon)) +
  geom_boxplot() +
  geom_point(aes(color = Sample), alpha = 0.7) +
  geom_line(aes(group = Sample, color = Sample), alpha = 0.7) +
  theme_bw() +
  theme(legend.position = "none")
  

kruskal.test(Shannon ~ method, data = stacked)
## BETA DIVERSITY ------------------------------------------------------------------------

robust_dist_methods <- vegdist(abund, method = "robust.aitchison")

permanova_methods <- adonis2(
  robust_dist_methods ~ method,
  data = stacked,
  permutations = 999
)

permanova_methods

library(vegan)
library(ggplot2)

# If you haven't already filtered out empty rows:
non_empty <- rowSums(abund) > 0
abund_clean <- abund[non_empty, , drop = FALSE]
stacked_clean <- stacked[non_empty, ]

robust_dist_methods <- vegdist(abund_clean, method = "robust.aitchison")

# Classical multidimensional scaling (PCoA)
ord <- cmdscale(robust_dist_methods, k = 2, eig = TRUE)

ord_df <- data.frame(
  MDS1   = ord$points[, 1],
  MDS2   = ord$points[, 2],
  method = stacked_clean$method,
  Sample = stacked_clean$Sample
)

table(ord_df$method)

ggplot(ord_df, aes(MDS1, MDS2, color = method)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(x = "PCoA1", y = "PCoA2", title = "Robust Aitchison PCoA by decontamination method")

ggplot(ord_df, aes(MDS1, MDS2, color = method, group = Sample)) +
  geom_point(alpha = 0.9
  )+
  geom_path(alpha = 0.3) +
  theme_minimal() +
  labs(x = "PCoA1", y = "PCoA2", title = "Sample trajectories across methods")

ggplot(ord_df, aes(MDS1, MDS2, color = method, group = Sample)) +
  geom_point(
    alpha = 0.9,
    position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_path(alpha = 0.3) +
  theme_minimal() +
  labs(x = "PCoA1", y = "PCoA2", title = "Sample trajectories across methods")



# SAVES THE DESIRED PHYLOSEQ OBJ (BEST DECONTAMINATION STRATEGY)

# merges with ps.raw (after removing shipment 5 samples from ps.raw)
ps.raw_no5 <- subset_samples(ps.raw, Shipment != 5)
ps.decontaminated <- merge_phyloseq(ps.raw_no5, ps.decontam_freandprev)

save(ps.decontaminated, file = "./Scripts/12s/rdata/ps.decontaminated.RData")
