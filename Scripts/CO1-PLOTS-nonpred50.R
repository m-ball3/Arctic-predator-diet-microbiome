# ------------------------------------------------------------------
# nonpred50 filtered plots
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
load("ps.CO1.nonpred50.Rdata")

ps.CO1 <- ps.CO1.nonpred50

## Merges same species
ps.CO1 = tax_glom(ps.CO1, "Species", NArm = FALSE)%>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Plots stacked bar plot of absolute abundance
plot_bar(ps.CO1, x="Specimen.ID", fill="Species")

# Calculates relative abundance of each species 
## Transforms NaN (0/0) to 0
psCO1.rel.nomams <- transform_sample_counts(ps.CO1, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(psCO1.rel.nomams))), arr.ind = TRUE)

# ------------------------------------------------------------------
# PLOTS BEFORE MAMMALIA FILTERING
# ------------------------------------------------------------------
# Creates bar plot of relative abundance
# Plots with WADE IDs - Species
sp.rel.plot <- plot_bar(psCO1.rel, fill="Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
sp.rel.plot

# Plots with ADFG IDs - Species
ADFG.sp <- plot_bar(psCO1.rel, x = "Specimen.ID", fill="Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
ADFG.sp

# Plots with WADE IDs - Genus
gen.rel.plot <- plot_bar(psCO1.rel, fill="Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
gen.rel.plot

# Plots with ADFG IDs - Genus
ADFG.gen <- plot_bar(psCO1.rel, x = "Specimen.ID", fill="Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
ADFG.gen

# Plots with WADE IDs - Family
fam.rel.plot <- plot_bar(psCO1.rel, fill="Family") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
fam.rel.plot

# Plots with ADFG IDs - Family
ADFG.fam <- plot_bar(psCO1.rel, x = "Specimen.ID", fill="Family") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )
ADFG.fam

# Facet wrapped by Sample_type - Species (WADE IDs)
faucet.samtype <- plot_bar(psCO1.rel, fill = "Species") +
  facet_wrap(~ Sample_type, ncol = 2, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 
faucet.samtype

# Facet wrapped by Sample_type - Species (ADFG IDs)
ADFG.faucet.samtype <- plot_bar(psCO1.rel, x = "Specimen.ID", fill = "Species") +
  facet_wrap(~ Sample_type, ncol = 2, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 
ADFG.faucet.samtype

# Facet wrapped by pred
faucet.pred <- plot_bar(psCO1.rel, fill = "Species") +
  facet_wrap(~ Predator, ncol = 2, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 
faucet.pred

# Facet wrapped by Sample_type - Species (ADFG IDs)
ADFG.faucet.pred <- plot_bar(psCO1.rel, x = "Specimen.ID", fill = "Species") +
  facet_wrap(~ Predator, ncol = 2, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "lines"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) 
ADFG.faucet.pred

#saves plots 
ggsave("Deliverables/CO1/nonpred50/WADE labels/CO1-speciesnonpred50.png", plot = sp.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/CO1/nonpred50/ADFG-CO1-speciesnonpred50.png", plot = ADFG.sp, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/CO1/nonpred50/WADE labels/CO1-genusnonpred50.png", plot = gen.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/CO1/nonpred50/ADFG-CO1-genusnonpred50.png", plot = ADFG.gen, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/CO1/nonpred50/WADE labels/CO1-familynonpred50.png", plot = fam.rel.plot, width = 16, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/CO1/nonpred50/ADFG-CO1-familynonpred50.png", plot = ADFG.fam, width = 16, height = 8, units = "in", dpi = 300)

ggsave("Deliverables/CO1/nonpred50/WADE labels/CO1-species-by-samtypenonpred50.png", plot = faucet.samtype, width = 30, height = 16, units = "in", dpi = 300)
ggsave("Deliverables/CO1/nonpred50/ADFG-CO1-species-by-samtypenonpred50.png", plot = ADFG.faucet.samtype, width = 30, height = 16, units = "in", dpi = 300)

ggsave("Deliverables/CO1/nonpred50/WADE labels/CO1-species-by-prednonpred50.png", plot = faucet.pred, width = 30, height = 16, units = "in", dpi = 300)
ggsave("Deliverables/CO1/nonpred50/ADFG-CO1-species-by-prednonpred50.png", plot = ADFG.faucet.pred, width = 30, height = 16, units = "in", dpi = 300)


# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(psCO1.rel.nomams))

# Changes NA.1 to it's corresponding ASV
taxa.names <- as.data.frame(tax_table(psCO1.rel.nomams)) %>% 
  rownames_to_column("ASV") %>% 
  mutate(Species = case_when(is.na(Species)~ASV,
                             TRUE~Species)) %>% 
  pull(Species)

colnames(otu.abs) <- taxa.names

tax_table <- as.data.frame(tax_table(psCO1.rel.nomams))

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
          "./Deliverables/CO1/nonpred50/ADFG_CO1_absolute_speciesxsamples-nonpred50.csv", 
          row.names = FALSE)

# Add Lab ID to otu.abs
write.csv(otu.prop %>% 
            rownames_to_column("LabID"), 
          "./Deliverables/CO1/nonpred50/ADFG_CO1_relative_speciesxsamples-nonpred50.csv", 
          row.names = FALSE)

write.csv(tax_table%>% 
            rownames_to_column("ASV"), 
          "./Deliverables/CO1/nonpred50/ADFG_CO1_tax_table-nonpred50.csv", row.names = FALSE)
