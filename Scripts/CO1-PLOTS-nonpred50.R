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
psCO1.rel <- transform_sample_counts(ps.CO1, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

#Checks for NaN's 
which(is.nan(as.matrix(otu_table(psCO1.rel))), arr.ind = TRUE)

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

