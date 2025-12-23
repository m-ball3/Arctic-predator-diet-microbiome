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

# Creates a vector of read counts
reads.arctic <- sample_sums(ps.12s.arctic)

# Gets the desired sample's read counts
reads.WADE.111 <- reads.arctic["WADE-003-111"]
reads.WADE.124 <- reads.arctic["WADE-003-124"]

reads.WADE.111 # 38402
reads.WADE.124 # 20768

# Transforms read counts to relative abundance of each species 
## Transforms NaN (0/0) to 0
ps12s.replicate.rel <- transform_sample_counts(ps.12s.replicates.arctic, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Create the stacked bar plot
plot_bar(ps12s.replicate.rel, x = "LabID", fill = "Species")

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
ps.12s.sbering <- prune_samples(sample_sums(ps.12s.sbering) >= 100, ps.12s.sbering)
sample_sums(ps.12s.sbering)
nsamples(ps.12s.sbering)

ps.12s.arctic <- prune_samples(sample_sums(ps.12s.arctic) >= 100, ps.12s.arctic)
sample_sums(ps.12s.arctic)
nsamples(ps.12s.arctic)

# ------------------------------------------------------------------
# EXPLORES SAMPLES LOST IN FILTERING FOR COOK INLET BELUGAS
# ------------------------------------------------------------------

# Keeps a copy before filtering
ps.12s.cook.before <- ps.12s.cook

# Applies filtering and keeeps a copy for after filtering
ps.12s.cook.after <- prune_samples(
  sample_sums(ps.12s.cook.before) >= 100,
  ps.12s.cook.before
)

# Breakdown of which samples are lost in filtering and their read counts


# Calculates the read counts before and after
cook.before.counts <- sample_sums(ps.12s.cook.before)
cook.after.counts  <- sample_sums(ps.12s.cook.after)

# Adds the sample names
cook.samples.before <- names(cook.before.counts)
cook.samples.after  <- names(cook.after.counts)

# Gets the samples that are lost iiin filtering 
cook.samples.lost <- setdiff(cook.samples.before, cook.samples.after)

# Creates a df to check what is lost in filtering\
# displays and gets mean and median
lost.df <- tibble(
  Sample = cook.samples.lost,
  Reads  = as.numeric(cook.before.counts[cook.samples.lost]))

lost.df
mean(lost.df$Reads)
median(lost.df$Reads)
range(lost.df$Reads)

# Plots a histogram of sample read count
p_lost_hist <- ggplot(lost.df, aes(x = Reads)) +
  geom_histogram(binwidth = 10, color = "black", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Cook – read counts of samples lost (before filter)",
    x = "Reads per sample",
    y = "Frequency"
  )
p_lost_hist

# Gets relative abundance for each
ps.12s.cook.before.rel <- transform_sample_counts(
  ps.12s.cook.before, function(x) x / sum(x))

ps.12s.cook.after.rel <- transform_sample_counts(
  ps.12s.cook.after, function(x) x / sum(x))

# Plots 
p_before_rel <- plot_bar(ps.12s.cook.before.rel, fill = "Species") +
  ggtitle("Cook – relative abundance BEFORE filter") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_after_rel <- plot_bar(ps.12s.cook.after.rel, fill = "Species") +
  ggtitle("Cook – relative abundance AFTER filter (≥100 reads)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# patchwork plot with one legend
(p_before_rel / p_after_rel +
    plot_layout(guides = "collect")) &
  theme(legend.position = "right")
# warning message related to removal of taxa with abundance of 0 


# ------------------------------------------------------------------
# CONTINUES WITH CLEANING PHYLOSEQ
# ------------------------------------------------------------------

## MERGE TO SPECIES HERE (TAX GLOM)
ps.12s.cook = tax_glom(ps.12s.cook.after, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

ps.12s.sbering = tax_glom(ps.12s.sbering, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

ps.12s.arctic = tax_glom(ps.12s.arctic, "Species", NArm = FALSE) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Filtering to remove taxa with less than 1% of reads assigned in at least 1 sample.
f1 <- filterfun_sample(function(x) x / sum(x) > 0.01)

lowcount.filt.cook <- genefilter_sample(ps.12s.cook.after, f1, A=1)
ps.12s.cook.filt <- prune_taxa(lowcount.filt.cook, ps.12s.cook.after)

lowcount.filt.sbering <- genefilter_sample(ps.12s.sbering, f1, A=1)
ps.12s.sbering.filt <- prune_taxa(lowcount.filt.sbering, ps.12s.sbering)

lowcount.filt.arctic <- genefilter_sample(ps.12s.arctic, f1, A=1)
ps.12s.arctic.filt <- prune_taxa(lowcount.filt.arctic, ps.12s.arctic)


# Explores ASV assignments
asv.cook.rows <- as.data.frame(tax_table(ps.12s.cook.filt))
asv.bering.rows <- as.data.frame(tax_table(ps.12s.sbering.filt))
asv.arctic.rows <- as.data.frame(tax_table(ps.12s.arctic.filt))

# makes rownames a column called rn
asv.cook.rows$rn <- rownames(asv.cook.rows)
asv.bering.rows$rn <- rownames(asv.bering.rows)
asv.bering.rows$rn <- rownames(asv.bering.rows)

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


# ------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------
# Creates bar plot of relative abundance

## SPECIES-LEVEL PLOTS ----------------------------------------------------

# Plots with WADE IDs
sp.cook.rel.plot <- plot_bar(ps12s.cook.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.cook.rel.plot

# Plots with ADFG IDs
sp.cook.rel.plot.ADFG.sp<- plot_bar(ps12s.cook.rel, x = "Specimen.ID", fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.cook.rel.plot.ADFG.sp

sp.sbering.rel.plot <- plot_bar(ps12s.sbering.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.sbering.rel.plot

# Plots with ADFG IDs
sp.sbering.rel.plot.ADFG.sp <- plot_bar(ps12s.sbering.rel, x = "Specimen.ID", fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.sbering.rel.plot.ADFG.sp

sp.arctic.rel.plot <- plot_bar(ps12s.arctic.rel, fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.arctic.rel.plot

# Plots with ADFG IDs
sp.arctic.rel.plot.ADFG.sp<- plot_bar(ps12s.arctic.rel, x = "Specimen.ID", fill="Species")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
sp.arctic.rel.plot.ADFG.sp

## GENUS-LEVEL PLOTS ----------------------------------------------------

# Cook Inlet
gen.cook.rel.plot <- plot_bar(ps12s.cook.rel, fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.cook.rel.plot

gen.cook.rel.plot.ADFG.gen <- plot_bar(ps12s.cook.rel, x = "Specimen.ID", fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.cook.rel.plot.ADFG.gen

# SE Bering
gen.sbering.rel.plot <- plot_bar(ps12s.sbering.rel, fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.sbering.rel.plot

gen.sbering.rel.plot.ADFG.gen <- plot_bar(ps12s.sbering.rel, x = "Specimen.ID", fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.sbering.rel.plot.ADFG.gen

# Arctic
gen.arctic.rel.plot <- plot_bar(ps12s.arctic.rel, fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.arctic.rel.plot

gen.arctic.rel.plot.ADFG.gen <- plot_bar(ps12s.arctic.rel, x = "Specimen.ID", fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gen.arctic.rel.plot.ADFG.gen


## FAMILY-LEVEL PLOTS ---------------------------------------------------

# Cook Inlet
fam.cook.rel.plot <- plot_bar(ps12s.cook.rel, fill = "Family") +
  theme_minimal() +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.cook.rel.plot

fam.cook.rel.plot.ADFG.fam <- plot_bar(ps12s.cook.rel, x = "Specimen.ID", fill = "Family") +
  theme_minimal() +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.cook.rel.plot.ADFG.fam

# SE Bering
fam.sbering.rel.plot <- plot_bar(ps12s.sbering.rel, fill = "Family") +
  theme_minimal() +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.sbering.rel.plot

fam.sbering.rel.plot.ADFG.fam <- plot_bar(ps12s.sbering.rel, x = "Specimen.ID", fill = "Family") +
  theme_minimal() +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.sbering.rel.plot.ADFG.fam

# Arctic
fam.arctic.rel.plot <- plot_bar(ps12s.arctic.rel, fill = "Family") +
  theme_minimal() +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.arctic.rel.plot

fam.arctic.rel.plot.ADFG.fam <- plot_bar(ps12s.arctic.rel, x = "Specimen.ID", fill = "Family") +
  theme_minimal() +
  labs(y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
fam.arctic.rel.plot.ADFG.fam


## FACETED BY PREDATOR (REGIONAL OBJECTS) ------------------------------

# WADE IDs
faucet.cook   <- plot_bar(ps12s.cook.rel,   x = "LabID",       fill = "Species") +
  facet_wrap(~ Predator, ncol = 3, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10))
  )
faucet.cook

faucet.sbering <- plot_bar(ps12s.sbering.rel, x = "LabID",     fill = "Species") +
  facet_wrap(~ Predator, ncol = 3, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10))
  )
faucet.sbering

faucet.arctic  <- plot_bar(ps12s.arctic.rel,  x = "LabID",     fill = "Species") +
  facet_wrap(~ Predator, ncol = 3, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10))
  )
faucet.arctic 

# ADFG IDs
faucet.cook.ADFG   <- plot_bar(ps12s.cook.rel,   x = "Specimen.ID", fill = "Species") +
  facet_wrap(~ Predator, ncol = 3, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10))
  )
faucet.cook.ADFG 

faucet.sbering.ADFG <- plot_bar(ps12s.sbering.rel, x = "Specimen.ID", fill = "Species") +
  facet_wrap(~ Predator, ncol = 3, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10))
  )
faucet.sbering.ADFG

faucet.arctic.ADFG  <- plot_bar(ps12s.arctic.rel,  x = "Specimen.ID", fill = "Species") +
  facet_wrap(~ Predator, ncol = 3, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0.5, "lines"),
    axis.title.x     = element_text(margin = margin(t = 10))
  )
faucet.arctic.ADFG

## SAVE PLOTS -----------------------------------------------------------

## Species-level, WADE IDs
ggsave("Deliverables/12S/regions/cook/WADE labels/12S-species-cook.png",
       plot = sp.cook.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/sbering/WADE labels/12S-species-sbering.png",
       plot = sp.sbering.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/arctic/WADE labels/12S-species-arctic.png",
       plot = sp.arctic.rel.plot, width = 16, height = 8, units = "in", dpi = 300)

## Species-level, ADFG IDs
ggsave("Deliverables/12S/regions/cook/ADFG-12S-species-cook.png",
       plot = sp.cook.rel.plot.ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/sbering/ADFG-12S-species-sbering.png",
       plot = sp.sbering.rel.plot.ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/arctic/ADFG-12S-species-arctic.png",
       plot = sp.arctic.rel.plot.ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)


## Genus-level, WADE IDs
ggsave("Deliverables/12S/regions/cook/WADE labels/12S-genus-cook.png",
       plot = gen.cook.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/sbering/WADE labels/12S-genus-sbering.png",
       plot = gen.sbering.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/arctic/WADE labels/12S-genus-arctic.png",
       plot = gen.arctic.rel.plot, width = 16, height = 8, units = "in", dpi = 300)

## Genus-level, ADFG IDs
ggsave("Deliverables/12S/regions/cook/ADFG-12S-genus-cook.png",
       plot = gen.cook.rel.plot.ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/sbering/ADFG-12S-genus-sbering.png",
       plot = gen.sbering.rel.plot.ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/arctic/ADFG-12S-genus-arctic.png",
       plot = gen.arctic.rel.plot.ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)


## Family-level, WADE IDs
ggsave("Deliverables/12S/regions/cook/WADE labels/12S-family-cook.png",
       plot = fam.cook.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/sbering/WADE labels/12S-family-sbering.png",
       plot = fam.sbering.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/arctic/WADE labels/12S-family-arctic.png",
       plot = fam.arctic.rel.plot, width = 16, height = 8, units = "in", dpi = 300)

## Family-level, ADFG IDs
ggsave("Deliverables/12S/regions/cook/ADFG-12S-family-cook.png",
       plot = fam.cook.rel.plot.ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/sbering/ADFG-12S-family-sbering.png",
       plot = fam.sbering.rel.plot.ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/12S/regions/arctic/ADFG-12S-family-arctic.png",
       plot = fam.arctic.rel.plot.ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)


## Faceted species-by-predator, WADE IDs
ggsave("Deliverables/12S/regions/cook/WADE labels/12S-species-by-pred-cook.png",
       plot = faucet.cook, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions/sbering/WADE labels/12S-species-by-pred-sbering.png",
       plot = faucet.sbering, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions/arctic/WADE labels/12S-species-by-pred-arctic.png",
       plot = faucet.arctic, width = 16, height = 8, units = "in", dpi = 300)

## Faceted species-by-predator, ADFG IDs
ggsave("Deliverables/12S/regions/cook/ADFG-12S-species-by-pred-cook.png",
       plot = faucet.cook.ADFG, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions//sbering/ADFG-12S-species-by-pred-sbering.png",
       plot = faucet.sbering.ADFG, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/12S/regions/arctic/ADFG-12S-species-by-pred-arctic.png",
       plot = faucet.arctic.ADFG, width = 16, height = 8, units = "in", dpi = 300)


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


## Cereates a function to convert absolute to relative by row
make_relative <- function(otu_abs) {
  otu_counts <- otu_abs[, -1]                     # removes Specimen.ID
  otu_prop   <- otu_counts / rowSums(otu_counts)  # calculates relative abundance
  otu_prop$Specimen.ID <- otu_abs$Specimen.ID     # adds Specimen.ID back in
  otu_prop <- otu_prop[, c(ncol(otu_prop), 1:(ncol(otu_prop)-1))]
  otu_prop[is.na(otu_prop)] <- 0                  # replaces NaNs with 0
  is.num <- sapply(otu_prop, is.numeric)
  otu_prop[is.num] <- lapply(otu_prop[is.num], round, 3)
  otu_prop
}

otu.cook.prop    <- make_relative(otu.cook.abs)
otu.sbering.prop <- make_relative(otu.sbering.abs)
otu.arctic.prop  <- make_relative(otu.arctic.abs)


## WRITE ABSOLUTE + RELATIVE TABLES AND TAX TABLES 

# Cook Inlet
write.csv(
  otu.cook.abs %>% rownames_to_column("LabID"),
  "./Deliverables/12S/regions/cook/ADFG_12s_cook_absolute_speciesxsamples.csv", row.names = FALSE)
write.csv(
  otu.cook.prop %>% rownames_to_column("LabID"),
  "./Deliverables/12S/regions/cook/ADFG_12s_cook_relative_speciesxsamples.csv", row.names = FALSE)
write.csv(
  cooktax_table %>% rownames_to_column("ASV"),
  "./Deliverables/12S/regions/cook/ADFG_12s_cook_tax_table.csv", row.names = FALSE)


# SE Bering
write.csv(
  otu.sbering.abs %>% rownames_to_column("LabID"),
  "./Deliverables/12S/regions/sbering/ADFG_12s_sbering_absolute_speciesxsamples.csv", row.names = FALSE)
write.csv(
  otu.sbering.prop %>% rownames_to_column("LabID"),
  "./Deliverables/12S/regions/sbering/ADFG_12s_sbering_relative_speciesxsamples.csv", row.names = FALSE)
write.csv(
  sberingtax_table %>% rownames_to_column("ASV"),
  "./Deliverables/12S/regions/sbering/ADFG_12s_sbering_tax_table.csv", row.names = FALSE)


# Arctic
write.csv(
  otu.arctic.abs %>% rownames_to_column("LabID"),
  "./Deliverables/12S/regions/arctic/ADFG_12s_arctic_absolute_speciesxsamples.csv", row.names = FALSE)
write.csv(
  otu.arctic.prop %>% rownames_to_column("LabID"),
  "./Deliverables/12S/regions/arctic/ADFG_12s_arctic_relative_speciesxsamples.csv", row.names = FALSE)
write.csv(
  arctictax_table %>% rownames_to_column("ASV"),
  "./Deliverables/12S/regions/arctic/ADFG_12s_arctic_tax_table.csv", row.names = FALSE)
