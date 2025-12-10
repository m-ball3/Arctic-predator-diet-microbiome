# ------------------------------------------------------------------
# 12S Samples by Location
# Mollie Ball
# 
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Sets up the Environment and Loads in data
# ------------------------------------------------------------------

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# 
# devtools::install_github("benjjneb/dada2", ref="v1.16", lib = .libPaths()[1])
# BiocManager::install("dada2", lib = .libPaths()[1], force = TRUE)
# BiocManager::install("S4Vectors")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(dplyr)
library(dada2)

# Loads in dada2 output
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output.Rdata")

# loads in regional DBs
cookinletDB <- "DADA2/Ref-DB/12S/12S_Cook-Inlet-DB.fasta"
cookinletDB.sp <- "DADA2/Ref-DB/12S/12S_Cook-Inlet-addspecies-DB.fasta"

sberingDB <- "DADA2/Ref-DB/12S/12S_S-Bering-DB.fasta"
sberingDB.sp <- "DADA2/Ref-DB/12S/12S_S-Bering-addspecies-DB.fasta"

arcticDB <- "DADA2/Ref-DB/12S/12S_Arctic-DB.fasta"
arcticDB.sp<- "DADA2/Ref-DB/12S/12S_Arctic-addspecies-DB.fasta"


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

cooktaxa <- assignTaxonomy(cook.seqtab , cookinletDB, tryRC = TRUE, minBoot = 95)
cooksp <- assignSpecies(cook.seqtab, cookinletDB.sp)%>%
  as.data.frame()%>%
  dplyr::rename(
    Genus.x = Genus, 
    Species.y = Species)%>%
  mutate(DB = "cookinletDB")%>%
  as.matrix()


sberingtaxa <- assignTaxonomy(sbering.seqtab, sberingDB, tryRC = TRUE, minBoot = 95)
sberingsp <- assignSpecies(sbering.seqtab, sberingDB.sp)%>%
  as.data.frame()%>%
  dplyr::rename(
    Genus.x = Genus, 
    Species.y = Species)%>%
  mutate(DB = "sberingDB")%>%
  as.matrix()

arctictaxa <- assignTaxonomy(arctic.seqtab, arcticDB, tryRC = TRUE, minBoot = 95)
arcticsp <- assignSpecies(arctic.seqtab, arcticDB.sp) %>%
  as.data.frame() %>%              
  dplyr::rename(
    Genus.x   = Genus,
    Species.y = Species) %>%
  mutate(DB = "arcticDB")%>%
  as.matrix()                      


# Combines taxa tables
taxam <- rbind(cooktaxa, sberingtaxa, arctictaxa)

# Combines species tables
genus.species <- rbind(cooksp, sberingsp, arcticsp)

taxa_na_fixed <- as.data.frame(taxam) %>% 
  rownames_to_column("ASV") %>%
  filter(is.na(Species)) %>%
  left_join(
    as.data.frame(genus.species) %>%
      rownames_to_column("ASV"),
    by = "ASV"
  ) %>%
  unite(col = addSpecies, Genus.x, Species.y, sep = " ") %>%
  ungroup() %>%
  mutate(addSpecies = case_when(
    addSpecies == "NA NA" ~ NA_character_,
    TRUE ~ addSpecies
  )) %>%
  mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%
  mutate(.grp = ifelse(is.na(addSpecies),
                       paste0("NA_grp_", row_number()),
                       addSpecies)) %>%
  group_by(.grp) %>%
  mutate(Class  = if (length(unique(Class))  > 1) NA else Class) %>%
  mutate(Order  = if (length(unique(Order))  > 1) NA else Order) %>%
  mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>%
  mutate(Genus  = if (length(unique(Genus))       > 1) NA else Genus) %>%
  ungroup() %>%
  select(-Species)%>%
  dplyr::rename(Species = addSpecies)

taxa <- bind_rows(
  taxa_na_fixed,
  as.data.frame(taxam) %>%
    rownames_to_column("ASV") %>%
    filter(!is.na(Species))
) %>%
  mutate(.grp = ifelse(is.na(Species),
                       paste0("NA_grp_", row_number()),
                       Species)) %>%
  group_by(.grp) %>%
  fill(Order, Family, Genus, .direction = "updown") %>%
  ungroup() %>%
  select(-.grp) %>%
  column_to_rownames("ASV") %>%
  as.matrix()

# Resaves output

save(seqtab.nochim, freq.nochim, track, taxa, file = "WADE003-arcticpred_dada2_QAQC_12SP1_output-regionalDB.Rdata")

