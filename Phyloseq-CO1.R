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
library(tibble)


# Loads dada2 output
#load("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")
load("DADA2/DADA2 Outputs/testDADA2_CO1_allseqs_012826.Rdata")

# ------------------------------------------------------------------
# FORMATS METADATASHEET FOR PHYLOSEQ OBJ
# ------------------------------------------------------------------
# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- sub("^((WADE-003-\\d+|WADE-003-\\d+-C|WADE-003-\\d+-UC))_.*", "\\1", rownames(seqtab.nochim))
rownames(seqtab.nochim) <- gsub("-16S_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv") %>%
  dplyr::rename(Predator = Species) %>% # renames species column to predator
  mutate(Predator = tolower(Predator))  # changes capitalization to all lowercase (fixes Beluga and beluga)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% 
      dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# Removes rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ]

# Sets row names to LabID
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

# ------------------------------------------------------------------
# CREATES PHYLOSEQ
# ------------------------------------------------------------------

# Creates master phyloseq object
ps.CO1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

nsamples(ps.CO1)

# Saves phyloseq obj
saveRDS(ps.CO1, "ps.CO1.raw")

# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES (NEW AND REPLACEMENTS)
# ------------------------------------------------------------------

# Identifies rows in Specimen.ID that appear more than once (replicates)
duplicated_ids <- samdf$Specimen.ID[duplicated(samdf$Specimen.ID) | duplicated(samdf$Specimen.ID, fromLast = TRUE)]
unique_dup_ids <- unique(duplicated_ids) 

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
ps.CO1.replicates <- subset_samples(ps.CO1, Specimen.ID %in% unique_dup_ids)
unique_dup_ids <- as.data.frame(unique_dup_ids)


## NEW AND REPLACEMENTS

# 08BS7 = WADE 094 and 086
# 09WWBS10 = WADE 096 and 085
# 2014BS18 = WADE 101-UC and 081
# 2021BDL-0723B = WADE 093 and 080
# 2023Beluga0721SA = WADE 099 and 077
# 2023Beluga0722SA = WADE 097 and 075
# 2023Beluga0722SF = 098 and 071
# DL22OTZ004 = WADE 092 and 040
# DL22SCM001 = WADE 091 and 041
# 2010RSW-05 = WADE 100-C, 100-UC, and 084
# 2014-02-RS = WADE 095 and 082

## Tech Rep

# PH22SH036-S = WADE 115 and 123
# PH23SH032-F = WADE 143 and 147

# Removes species assignments less than 100 reads for readability
ps.CO1.replicates <- prune_taxa(taxa_sums(ps.CO1.replicates) > 0, ps.CO1.replicates)

# Creates a vector of read counts
reads <- sample_sums(ps.CO1)

# Gets the desired sample's read counts
# NEW AND REPLACEMENTS -----------------------------------------

# 08BS7 = WADE 094 and 086
reads.WADE.094 <- reads["WADE-003-094"] # 4592
reads.WADE.086 <- reads["WADE-003-086"] # 180

# 09WWBS10 = WADE 096 and 085
reads.WADE.096 <- reads["WADE-003-096"] # 31866
reads.WADE.085 <- reads["WADE-003-085"] # 11007

# 2014BS18 = WADE 101-UC and 081
reads.WADE.101UC <- reads["WADE-003-101-UC"] # 1312
reads.WADE.081   <- reads["WADE-003-081"]    # 7017

# 2021BDL-0723B = WADE 093 and 080
reads.WADE.093 <- reads["WADE-003-093"] # 10
reads.WADE.080 <- reads["WADE-003-080"] # 2955

# 2023Beluga0721SA = WADE 099 and 077
reads.WADE.099 <- reads["WADE-003-099"] # 9354
reads.WADE.077 <- reads["WADE-003-077"] # 1414

# 2023Beluga0722SA = WADE 097 and 075
reads.WADE.097 <- reads["WADE-003-097"] # 8907
reads.WADE.075 <- reads["WADE-003-075"] # 3247

# 2023Beluga0722SF = 098 and 071
reads.WADE.098 <- reads["WADE-003-098"] # 9882
reads.WADE.071 <- reads["WADE-003-071"] # 3779

# DL22OTZ004 = WADE 092 and 040
reads.WADE.092 <- reads["WADE-003-092"] # 13086
reads.WADE.040 <- reads["WADE-003-040"] # 147

# DL22SCM001 = WADE 091 and 041
reads.WADE.091 <- reads["WADE-003-091"] # 17625
reads.WADE.041 <- reads["WADE-003-041"] # 996

# 2010RSW-05 = WADE 100-C, 100-UC, and 084
reads.WADE.100C  <- reads["WADE-003-100-C"] # 3461
reads.WADE.100UC <- reads["WADE-003-100-UC"]# 30313
reads.WADE.084   <- reads["WADE-003-084"]   # 127

# 2014-02-RS = WADE 095 and 082
reads.WADE.095 <- reads["WADE-003-095"] # 21626
reads.WADE.082 <- reads["WADE-003-082"] # 28363


# TECH REP ------------------------------------------------------

# PH22SH036-S = WADE 115 and 123
reads.WADE.115 <- reads["WADE-003-115"] # 5096
reads.WADE.123 <- reads["WADE-003-123"] # 1685

# PH23SH032-F = WADE 143 and 147
reads.WADE.143 <- reads["WADE-003-143"] # 23930
reads.WADE.147 <- reads["WADE-003-147"] # 7124

# Drop Mammalia, keep all other classes (including NA)
ps.CO1.replicates <- subset_taxa(
  ps.CO1.replicates,
  Class != "Mammalia" | is.na(Class)
)

nsamples(ps.CO1.replicates)

# Creates a vector of read counts
reads <- sample_sums(ps.CO1.replicates)


# Create the stacked bar plot (absolute)
plot_bar(ps.CO1.replicates, x = "LabID", fill = "Genus") +
  facet_wrap(~ Specimen.ID, ncol = 13, scales = "free_x", strip.position = "top")

# Transforms read counts to relative abundance of each species 
## Transforms NaN (0/0) to 0
psCO1.replicate.rel <- transform_sample_counts(ps.CO1.replicates, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Create the stacked bar plot (relative)
plot_bar(psCO1.replicate.rel, x = "LabID", fill = "Genus")+
  facet_wrap(~ Specimen.ID, ncol = 13, scales = "free_x", strip.position = "top")

# Ensures samples removed in filtering are removed from samdf
replicate_to_remove <- c(
  "WADE-003-086",
  "WADE-003-085",
  "WADE-003-101-UC",
  "WADE-003-093",
  "WADE-003-077",
  "WADE-003-075",
  "WADE-003-071",
  "WADE-003-040",
  "WADE-003-041",
  "WADE-003-100-C",
  "WADE-003-084",
  "WADE-003-095",
  "WADE-003-123",
  "WADE-003-147"
)

samdf <- samdf[!rownames(samdf) %in% replicate_to_remove, ]

# Recreates phyloseq object without unwanted replicate
ps.CO1 <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
  sample_data(samdf), 
  tax_table(compiled_taxa)
)%>% 
  subset_samples(LabID != replicate_to_remove)
sample_names(ps.CO1)
nsamples(ps.CO1)

# ------------------------------------------------------------------
# CLEANS PHYLOSEQ WITH REMOVAL TRACKING
# ------------------------------------------------------------------

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.CO1))
names(dna) <- taxa_names(ps.CO1)
ps.raw <- merge_phyloseq(ps.CO1, dna)
taxa_names(ps.CO1) <- paste0("ASV", seq(ntaxa(ps.CO1)))

# Initialize removal tracking dataframe
removal_log <- data.frame(
  SampleID = character(),
  Step = character(),
  ReadCount = numeric(),
  stringsAsFactors = FALSE
)

# Step 2: Remove samples with total abundance < 100
samples_before <- sample_names(ps.CO1)
reads_before <- sample_sums(ps.CO1)
ps.CO1 <- prune_samples(sample_sums(ps.CO1) >= 100, ps.CO1)
low_abund_samples <- setdiff(samples_before, sample_names(ps.CO1))

# Log low abundance removals
if(length(low_abund_samples) > 0) {
  removal_log <- rbind(removal_log, 
                       data.frame(
                         SampleID = low_abund_samples,
                         Step = "Low abundance (<100 reads)",
                         ReadCount = reads_before[low_abund_samples]
                       )
  )
}

# Step 3: Prevalence filtering (taxa only, no samples removed)
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)
ps.CO1 <- prune_taxa(genefilter_sample(ps.CO1, f1, A=1), ps.CO1)

# Show removal summary
print("=== FILTERING REMOVAL SUMMARY ===")
print(removal_log)
print(paste("Total samples removed:", nrow(removal_log)))
print(paste("Samples remaining:", nsamples(ps.CO1)))

# # Save removal log
# write.csv(removal_log, "filtering_removal_log.csv", row.names = FALSE)

# Sync samdf
row_to_remove <- removal_log$SampleID
samdf <- samdf[!rownames(samdf) %in% row_to_remove, ]

# Filter out Mammalia (samples unaffected)
ps.CO1 <- subset_taxa(ps.CO1, class != "Mammalia" | is.na(class))
nsamples(ps.CO1)

# # Saves phyloseq obj
# saveRDS(ps.CO1, "ps.CO1")

## Merges same species
ps.CO1 = tax_glom(ps.CO1, "genus", NArm = FALSE)%>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Plots stacked bar plot of absolute abundance
plot_bar(ps.CO1, x="Specimen.ID", fill="genus")

# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
ps16s.rel <- transform_sample_counts(ps.CO1, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(ps16s.rel))), arr.ind = TRUE)

# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

# Plots with WADE IDs
sp.rel.plot <- plot_bar(ps16s.rel, fill="Species.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.rel.plot

# Plots with ADFG IDs
ADFG.sp<- plot_bar(ps16s.rel, x = "Specimen.ID", fill="Species.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.sp

gen.rel.plot <- plot_bar(ps16s.rel, fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.rel.plot

ADFG.gen <- plot_bar(ps16s.rel, x = 'Specimen.ID', fill="Genus.y")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.gen

fam.rel.plot <- plot_bar(ps16s.rel, fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.rel.plot

ADFG.fam <- plot_bar(ps16s.rel, "Specimen.ID", fill="Family")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ADFG.fam

# Facet wrapped by predator species
faucet <- plot_bar(ps16s.rel, fill = "Species.y") +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 

ADFG.faucet <- plot_bar(ps16s.rel, x ="Specimen.ID", fill = "Species.y") +
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
ggsave("Deliverables/16S/WADE labels/16S-species.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-species.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/WADE labels/16S-genus.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-genus.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/WADE labels/16S-family.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-family.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/16S/WADE labels/16S-species-by-pred.111125.png", plot = faucet, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/16S/ADFG-16S-species-by-pred.111125.png", plot = ADFG.faucet, width = 16, height = 8, units = "in", dpi = 300)

# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.CO1))

# Changes NA.1 to it's corresponding ASV
taxa.names <- as.data.frame(tax_table(ps.CO1)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species.y = case_when(is.na(Species.y)~ASV,
                               TRUE~Species.y)) %>% 
  pull(Species.y)

colnames(otu.abs) <- taxa.names

tax_table <- as.data.frame(tax_table(ps.CO1))

## Adds ADFG Sample ID as a column
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

# Add Lab ID to otu.abs
write.csv(otu.abs %>% 
            rownames_to_column("LabID"), 
          "./Deliverables/16S/ADFG_16s_absolute_speciesxsamples.csv", 
          row.names = FALSE)

# Add Lab ID to otu.abs
write.csv(otu.prop %>% 
            rownames_to_column("LabID"), 
          "./Deliverables/16S/ADFG_16s_relative_speciesxsamples.csv", 
          row.names = FALSE)

write.csv(tax_table%>% 
            rownames_to_column("ASV"), 
          "./Deliverables/16S/ADFG_16s_tax_table.csv", row.names = FALSE)


