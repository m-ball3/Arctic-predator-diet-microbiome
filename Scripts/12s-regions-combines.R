# ------------------------------------------------------------------
# TOWARDS COMBINING REGIONAL 12S PHYLOSEQ OBJS
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

load("ps.12s.regions.filt.Rdata")

# ------------------------------------------------------------------
# Deals with ASV assignments
# ------------------------------------------------------------------

# Adds marker to ASV to indicate DB
taxa_names(ps.12s.cook.filt)    <- paste0(taxa_names(ps.12s.cook.filt), ".CI")
taxa_names(ps.12s.sbering.filt) <- paste0(taxa_names(ps.12s.sbering.filt), ".SB")
taxa_names(ps.12s.arctic.filt)  <- paste0(taxa_names(ps.12s.arctic.filt), ".A")

# Checks
asv.cook.rows <- as.data.frame(tax_table(ps.12s.cook.filt))
asv.bering.rows <- as.data.frame(tax_table(ps.12s.sbering.filt))
asv.arctic.rows <- as.data.frame(tax_table(ps.12s.arctic.filt))


