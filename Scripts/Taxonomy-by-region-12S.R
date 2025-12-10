# ------------------------------------------------------------------
# 12S Samples by Location
# Mollie Ball
# 
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Sets up the Environment and Loads in data
# ------------------------------------------------------------------

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

devtools::install_github("benjjneb/dada2", ref="v1.16", lib = .libPaths()[1])
BiocManager::install("dada2", lib = .libPaths()[1], force = TRUE)
BiocManager::install("S4Vectors")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(dplyr)
library(dada2)

# Loads in dada2 output
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output.Rdata")

# loads in regional DBs
cookinletDB <- ".DADA2/Ref-DB/12S/12S_Cook-Inlet-DB.fasta"
cookinletDB.sp <- ".DADA2/Ref-DB/12S/12S_Cook-Inlet-addspecies-DB.fasta"

sberingDB <- ".DADA2/Ref-DB/12S/12S_S-Bering-DB.fasta"
sberingDB.sp <- ".DADA2/Ref-DB/12S/12S_S-Bering-addspecies-DB.fasta"

arcticDB <- ".DADA2/Ref-DB/12S/12S_Arctic-DB.fasta"
arcticDB.sp<- ".DADA2/Ref-DB/12S/12S_Arctic-addspecies-DB.fasta"


# ------------------------------------------------------------------
# FORMATS METADATASHEET FOR PHYLOSEQ OBJ
# ------------------------------------------------------------------

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-MFU_S\\d+", "", rownames(seqtab.nochim))

# Gets sample metadata; filters out NA's (shipment 1)
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!is.na(LabID))

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator" & makes all lowercase
samdf <- dplyr::rename(samdf, Predator = Species)
samdf$Predator <- tolower(samdf$Predator)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# Remove rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ]

# Then set row names to LabID
rownames(samdf) <- samdf$LabID

# Only keeps rows that appear in both metadata and seq.tab 
## AKA only samples that made it through all steps 
common_ids <- intersect(rownames(samdf), rownames(seqtab.nochim))
samdf <- samdf[common_ids, ]
seqtab.nochim <- seqtab.nochim[common_ids, ]

# Checks for duplicates and identical sample rownames in both

### should return FALSE
any(duplicated(rownames(samdf)))
any(duplicated(rownames(seqtab.nochim)))

### should return TRUE
all(rownames(samdf) %in% rownames(seqtab.nochim))
all(rownames(seqtab.nochim) %in% rownames(samdf))

# Checks for samples in metadata but not in OTU table
setdiff(rownames(samdf), rownames(seqtab.nochim))

# Checks for samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf))

# Sanity check: row names are the same
rownames(samdf)
rownames(seqtab.nochim)


# ------------------------------------------------------------------
# FORMATS FOR REGION-SPECIFIC ASSIGNMENT
# ------------------------------------------------------------------

# adds row specifying DB
samdf <- samdf %>%
  mutate(
    DB = case_when(
      Location == "Cook Inlet" ~ "cookinletDB",
      Location %in% c("Hooper Bay", "Scammon Bay") ~ "sberingDB",
      TRUE ~ "arcticDB"
    )
  )

# Divides seqtab.nochim by lab ID into regions for Assign Taxonomy and Species
cook.ids   <- samdf$LabID[samdf$DB == "cookinletDB"]
sbering.ids <- samdf$LabID[samdf$DB == "sberingDB"]
arctic.ids <- samdf$LabID[samdf$DB == "arcticDB"]

cook.seqtab   <- seqtab.nochim[rownames(seqtab.nochim) %in% cook.ids, ]
sbering.seqtab <- seqtab.nochim[rownames(seqtab.nochim) %in% sbering.ids, ]
arctic.seqtab  <- seqtab.nochim[rownames(seqtab.nochim) %in% arctic.ids, ]

  
# ------------------------------------------------------------------
# SPECIES ASSIGNMENTS BY REGION
# ------------------------------------------------------------------

# Assigns Taxonomy and Species

cooktaxa <- assignTaxonomy(seqtab.nochimcook, cookinletDB, tryRC = TRUE, minBoot = 95)
cooksp <- assignSpecies(seqtab.nochim, cookinletDB.sp)

sberingtaxa <- assignTaxonomy(seqtab.nochim, sberingDB, tryRC = TRUE, minBoot = 95)
sberingsp <- assignSpecies(seqtab.nochim, sberingDB.sp)

arctictaxa <- assignTaxonomy(seqtab.nochim, arcticDB, tryRC = TRUE, minBoot = 95)
arcticsp <-assignSpecies(seqtab.nochim, arcticDB.sp)

# Adds species to tax table
cook <- addSpecies(cooktaxa, cooksp, verbose=TRUE)

sbering <- addSpecies(sberingtaxa, sberingsp, verbose=TRUE)

arctic <- addSpecies(arctictaxa, arcticsp, verbose=TRUE)

# Combines all three tax tables (same columns, stacked by rows)
taxa <- rbind(cook, sbering, arctic)
  
# Resaves output

save(seqtab.nochim, freq.nochim, track, taxa, file = "WADE003-arcticpred_dada2_QAQC_12SP1_output-regionalDB.Rdata")

