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
load("C:/Users/Intern/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP2_output2.Rdata")

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-MFU_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata
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
samdf <- samdf[!is.na(samdf$LabID), ]

# Then set row names to LabID
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

# Creates stacked bar plot 
## proportion of each species 
### samples on x
ps12s.prop <- transform_sample_counts(ps.12s, function(otu) otu/sum(otu))
ps12s.bar <- transform_sample_counts(ps.12s, function(OTU) OTU/sum(OTU))
ps12s.bar <- prune_taxa(taxa_names(ps12s.prop), ps12s.bar)

plot_bar(ps12s.bar, x="LabID", fill="Species") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# FIlters NAs
ps12s.bar.no.na <- subset_taxa(ps12s.bar, !is.na(Species))

plot_bar(ps12s.bar.no.na, x = "LabID", fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# FIlters beluga whale as prey species
ps12s.bar.no.beluga <- subset_taxa(ps12s.bar.no.na, Species != "Delphinapterus leucas")

plot_bar(ps12s.bar.no.beluga, x = "LabID", fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Facet wrapped by predator species
plot_bar(ps12s.bar.no.beluga, x="LabID", fill="Species") +
  facet_wrap(~ Predator, ncol = 1, strip.position = "right") +
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
sample_data(ps12s.bar.no.beluga)$Stomach.Goo <- factor(
  sample_data(ps12s.bar.no.beluga)$Stomach.Goo,
  levels = c("No", "Yes"),
  labels = c("Fecal", "Stomach")
)


plot_bar(ps12s.bar.no.beluga, x = "LabID", fill = "Species") +
  facet_wrap(~ Stomach.Goo, ncol = 1, strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  guides(fill = guide_legend(title = "Species"))

### export csv for ampbias correction
# readcount.table <- as.data.frame(otu_table(ps.raw))
# taxon.table <- as.data.frame(tax_table(ps.raw)) %>% rownames_to_column(var = "ASV")
# metadata.table <- samdf %>% rownames_to_column(var = "Sample")
# reference.table <- as.data.frame(refseq(ps.raw)) %>% 
#   rownames_to_column("ASVname")
