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
load("C:/Users/Intern/Arctic-predator-diet-microbiome/DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP2_output.Rdata")

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-16S_S\\d+", "", rownames(seqtab.nochim))

# # Gets sample metadata
# labdf <- read.csv("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/metadata/ADFG_dDNA_labwork_metadata.csv")%>%
#   filter(!is.na(LabID))
# 
# samdf <- read.csv("C:/Users/MBall/OneDrive/文档/WADE LAB/Arctic-predator-diet-microbiome/metadata/ADFG_dDNA_sample_metadata.csv")

# Gets sample metadata
labdf <- read.csv("C:/Users/Intern/Arctic-predator-diet-microbiome/metadata/ADFG_dDNA_labwork_metadata.csv") %>%
   filter(!is.na(LabID))

samdf <- read.csv("C:/Users/Intern/Arctic-predator-diet-microbiome/metadata/ADFG_dDNA_sample_metadata.csv")


# 1. Create mapping table with BOTH Specimen.ID AND Repeat.or.New.Specimen
map_unique <- labdf %>%
  filter(!is.na(LabID)) %>%
  distinct(Specimen.ID, Repeat.or.New.Specimen., LabID)

# 2. Join to samdf using BOTH columns
samdf <- samdf %>%
  left_join(map_unique, by = c("Specimen.ID", "Repeat.or.New.Specimen.")) %>%
  filter(!is.na(LabID)) %>%
  column_to_rownames("LabID")

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)


# Creates master phyloseq object
ps.16s <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.16s))
names(dna) <- taxa_names(ps.16s)
ps.raw <- merge_phyloseq(ps.16s, dna)
taxa_names(ps.16s) <- paste0("ASV", seq(ntaxa(ps.16s)))

nsamples(ps.16s)

# Creates stacked bar plot 
## proportion of each species 
  ### samples on x
ps16s.prop <- transform_sample_counts(ps.16s, function(otu) otu/sum(otu))
ps16s.bar <- transform_sample_counts(ps.16s, function(OTU) OTU/sum(OTU))
ps16s.bar <- prune_taxa(taxa_names(ps16s.prop), ps16s.bar)

plot_bar(ps16s.bar, x="Specimen.ID", fill="Species") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

plot_bar(ps16s.bar, x="Specimen.ID", fill="Species") +
  facet_wrap("Specimen.ID") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))








### export csv for ampbias correction
readcount.table <- as.data.frame(otu_table(ps.raw))
taxon.table <- as.data.frame(tax_table(ps.raw)) %>% rownames_to_column(var = "ASV")
metadata.table <- samdf %>% rownames_to_column(var = "Sample")
reference.table <- as.data.frame(refseq(ps.raw)) %>% 
  rownames_to_column("ASVname")
