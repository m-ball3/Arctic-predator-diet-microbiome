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
## TEMPORARILY RENAMES EB22PH005-S TO EB23PH005-S -----------------------------------------------------------

labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!is.na(LabID))%>%
  mutate(Specimen.ID = ifelse(Specimen.ID == "EB22PH005-S", "EB23PH005-S", Specimen.ID))

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv") %>%
  dplyr::rename(Predator = Species) %>%
  mutate(Predator = tolower(Predator)) %>%
  left_join( # Creates a column corresponding ADFG sample IDs with WADE sample IDs
    labdf %>% dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  ) %>%
  filter(!is.na(LabID)) %>% # Removes rows where LabID is NA (because shipment 1 was bad & thus not extracted)
  mutate(
    LabID = gsub("-C$", "", LabID)   # Removes "-C" (all 12s samples are cleaned)
  )

#  Sets row names to LabID
rownames(samdf) <- samdf$LabID

# INVESTIGATES WHY WE ADD 10 OBS FROM CREATING A CORRESPONDING ADFG SAMPLE ID WITH WADE ID

# Check for duplicate keys in labdf
labdf_dupes <- labdf %>%
  dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID) %>%
  group_by(Specimen.ID, Repeat.or.New.Specimen.) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

print(labdf_dupes)  # Shows which keys have multiples
nrow(labdf_dupes)   # Number of duplicated keys (likely 10)

labdf_dupes_with_labid <- labdf %>%
  dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID) %>%
  inner_join(labdf_dupes, by = c("Specimen.ID", "Repeat.or.New.Specimen.")) %>%
  arrange(Specimen.ID, Repeat.or.New.Specimen.)

print(labdf_dupes_with_labid)

# All samples after DADA2, before metadata filtering
seq_samples <- rownames(seqtab.nochim)      # should be length 79

# All samples with metadata after LabID join and NA removal
meta_samples <- rownames(samdf)             # before intersect(), should be > 74

# Only keeps rows that appear in both metadata and seq.tab 
## AKA only samples that made it through all steps 
common_ids <- intersect(rownames(samdf), rownames(seqtab.nochim)) # CURRENTLY 78; resolve sample 122 field ID mismatch
samdf <- samdf[common_ids, ]
seqtab.nochim <- seqtab.nochim[common_ids, ]

# Samples that disappear after intersect()
dropped_from_seq <- setdiff(seq_samples, meta_samples)
dropped_from_meta <- setdiff(meta_samples, seq_samples)

dropped_from_seq # 122
dropped_from_meta # make sure none of these are removed in error

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
    sample_DB = case_when(
      Location == "Cook Inlet" ~ "cookinletDB",
      Location %in% c("Hooper Bay", "Scammon Bay") ~ "sberingDB",
      Predator == "beluga whale" & Location == "Nome" ~ "sberingDB",
      TRUE ~ "arcticDB"
    )
  )

# Divides seqtab.nochim by lab ID into regions for Assign Taxonomy and Species
rowcount <- rownames(seqtab.nochim) # 78 as of 12/23/25

cook.ids   <- samdf$LabID[samdf$sample_DB == "cookinletDB"]
sbering.ids <- samdf$LabID[samdf$sample_DB == "sberingDB"]
arctic.ids <- samdf$LabID[samdf$sample_DB == "arcticDB"]

cook.seqtab   <- seqtab.nochim[rownames(seqtab.nochim) %in% cook.ids, ]
sbering.seqtab <- seqtab.nochim[rownames(seqtab.nochim) %in% sbering.ids, ]
arctic.seqtab  <- seqtab.nochim[rownames(seqtab.nochim) %in% arctic.ids, ]

rowcount.cook <- rownames(cook.seqtab) # 19
rowcount.sbering <- rownames(sbering.seqtab) # 19
rowcount.arctic <- rownames(arctic.seqtab) # 40
# total = 78 as of 12/23/25

# ------------------------------------------------------------------
# SPECIES ASSIGNMENTS BY REGION
# ------------------------------------------------------------------

# Assigns Taxonomy and Species

cooktaxa <- assignTaxonomy(cook.seqtab , cookinletDB, tryRC = TRUE, minBoot = 95)%>%
  as.data.frame()%>%
  as.matrix()

cooksp <- assignSpecies(cook.seqtab, cookinletDB.sp)%>%
  as.data.frame()%>%
  dplyr::rename(
    Genus.x = Genus, 
    Species.y = Species)%>%
  as.matrix()


sberingtaxa <- assignTaxonomy(sbering.seqtab, sberingDB, tryRC = TRUE, minBoot = 95)%>%
  as.data.frame()%>%
  as.matrix()
  
sberingsp <- assignSpecies(sbering.seqtab, sberingDB.sp)%>%
  as.data.frame()%>%
  dplyr::rename(
    Genus.x = Genus, 
    Species.y = Species)%>%
  as.matrix()

arctictaxa <- assignTaxonomy(arctic.seqtab, arcticDB, tryRC = TRUE, minBoot = 95)%>%
  as.data.frame() %>%
  as.matrix() 
  
arcticsp <- assignSpecies(arctic.seqtab, arcticDB.sp) %>%
  as.data.frame() %>%              
  dplyr::rename(
    Genus.x   = Genus,
    Species.y = Species) %>%
  as.matrix()                      

# ------------------------------------------------------------------
# ADDRESSES NA AND COLUMN NON-AGREEMENT ISSUES BETWEEN ASSIGNTAXONOMY AND ADDSPECIES
# ------------------------------------------------------------------

# COOK INLET DB
taxa.cook <- as.data.frame(cooktaxa) %>% 
  rownames_to_column("ASV") %>%  
  filter(is.na(Species)) %>%  
  left_join(as.data.frame(cooksp) %>% 
              rownames_to_column("ASV"),  by = "ASV") %>%  
  unite(col = addSpecies, Genus.x, Species.y, sep = " ") %>%  
  ungroup() %>%  
  mutate(addSpecies = case_when(addSpecies == "NA NA"~NA, TRUE~addSpecies)) %>% 
  mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%  
  mutate(
    Genus_from_sp = case_when(
      grepl(" ", addSpecies) ~ word(addSpecies, 1),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Genus = ifelse(is.na(Genus), Genus_from_sp, Genus)) %>%  # Backfill Genus NAs
  mutate(.grp = ifelse(is.na(addSpecies),  
                       paste0("NA_grp_", 
                              row_number()),  
                       addSpecies)) %>%  
  group_by(.grp) %>%  
  mutate(Class = if (length(unique(Class)) > 1) NA else Class) %>% 
  mutate(Order = if (length(unique(Order)) > 1) NA else Order) %>% 
  mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>% 
  mutate(Genus.x = if (length(unique(Genus)) > 1) NA else Genus) %>% 
  ungroup() %>%  
  select(-Species, -Genus.x, -Genus_from_sp) %>% 
  dplyr::rename("Species" = addSpecies) %>%  
  bind_rows(as.data.frame(cooktaxa) %>% 
              rownames_to_column("ASV") %>%  
              filter(!is.na(Species))) %>%  
  mutate(.grp = ifelse(is.na(Species),  
                       paste0("NA_grp_", 
                              row_number()),  
                       Species)) %>%  
  group_by(.grp) %>%  
  fill(Order, Family, Genus, .direction = "updown") %>% 
  ungroup() %>%  
  select(-.grp) %>%   
  mutate(DB = "cookinletDB") %>%
  relocate(DB, .before = Kingdom) %>%
  column_to_rownames("ASV") %>%  
  mutate(Species = ifelse(is.na(Species) & !is.na(Genus),
                          paste0(Genus, " spp."),
                          Species)) %>%
  as.matrix()

# S BERING DB
taxa.sbering <- as.data.frame(sberingtaxa) %>% 
  rownames_to_column("ASV") %>%  
  filter(is.na(Species)) %>%  
  left_join(as.data.frame(sberingsp) %>% 
              rownames_to_column("ASV"),  by = "ASV") %>%  
  unite(col = addSpecies, Genus.x, Species.y, sep = " ") %>%  
  ungroup() %>%  
  mutate(addSpecies = case_when(addSpecies == "NA NA"~NA, TRUE~addSpecies)) %>% 
  mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%  
  mutate(
    Genus_from_sp = case_when(
      grepl(" ", addSpecies) ~ word(addSpecies, 1),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Genus = ifelse(is.na(Genus), Genus_from_sp, Genus)) %>%  # Backfill Genus NAs
  mutate(.grp = ifelse(is.na(addSpecies),  
                       paste0("NA_grp_", 
                              row_number()),  
                       addSpecies)) %>%  
  group_by(.grp) %>%  
  mutate(Class = if (length(unique(Class)) > 1) NA else Class) %>% 
  mutate(Order = if (length(unique(Order)) > 1) NA else Order) %>% 
  mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>% 
  mutate(Genus.x = if (length(unique(Genus)) > 1) NA else Genus) %>% 
  ungroup() %>%  
  select(-Species, -Genus.x, -Genus_from_sp) %>% 
  dplyr::rename("Species" = addSpecies) %>%  
  bind_rows(as.data.frame(sberingtaxa) %>% 
              rownames_to_column("ASV") %>%  
              filter(!is.na(Species))) %>%  
  mutate(.grp = ifelse(is.na(Species),  
                       paste0("NA_grp_", 
                              row_number()),  
                       Species)) %>%  
  group_by(.grp) %>%  
  fill(Order, Family, Genus, .direction = "updown") %>% 
  ungroup() %>%  
  select(-.grp) %>% 
  mutate(DB = "sberingDB") %>%
  relocate(DB, .before = Kingdom) %>%
  column_to_rownames("ASV") %>%  
  mutate(Species = ifelse(is.na(Species) & !is.na(Genus),
                          paste0(Genus, " spp."),
                          Species)) %>%
  as.matrix()

# ARCTIC DB
taxa.arctic <- as.data.frame(arctictaxa) %>% 
  rownames_to_column("ASV") %>%  
  filter(is.na(Species)) %>%  
  left_join(as.data.frame(arcticsp) %>% 
              rownames_to_column("ASV"),  by = "ASV") %>%  
  unite(col = addSpecies, Genus.x, Species.y, sep = " ") %>%  
  ungroup() %>%  
  mutate(addSpecies = case_when(addSpecies == "NA NA"~NA, TRUE~addSpecies)) %>% 
  mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%  
  mutate(
    Genus_from_sp = case_when(
      grepl(" ", addSpecies) ~ word(addSpecies, 1),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(Genus = ifelse(is.na(Genus), Genus_from_sp, Genus)) %>%  # Backfill Genus NAs
  mutate(.grp = ifelse(is.na(addSpecies),  
                       paste0("NA_grp_", 
                              row_number()),  
                       addSpecies)) %>%  
  group_by(.grp) %>%  
  mutate(Class = if (length(unique(Class)) > 1) NA else Class) %>% 
  mutate(Order = if (length(unique(Order)) > 1) NA else Order) %>% 
  mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>% 
  mutate(Genus.x = if (length(unique(Genus)) > 1) NA else Genus) %>% 
  ungroup() %>%  
  select(-Species, -Genus.x, -Genus_from_sp) %>% 
  dplyr::rename("Species" = addSpecies) %>%  
  bind_rows(as.data.frame(arctictaxa) %>% 
              rownames_to_column("ASV") %>%  
              filter(!is.na(Species))) %>%  
  mutate(.grp = ifelse(is.na(Species),  
                       paste0("NA_grp_", 
                              row_number()),  
                       Species)) %>%  
  group_by(.grp) %>%  
  fill(Order, Family, Genus, .direction = "updown") %>% 
  ungroup() %>%  
  select(-.grp) %>%  
  mutate(DB = "arcticDB") %>%
  relocate(DB, .before = Kingdom) %>%
  column_to_rownames("ASV") %>%  
  mutate(Species = ifelse(is.na(Species) & !is.na(Genus),
                          paste0(Genus, " spp."),
                          Species)) %>%
  as.matrix()

# Resaves output

save(samdf, cook.seqtab, freq.nochim, track, taxa.cook, file = "DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-COOKINLET.Rdata")
save(samdf, sbering.seqtab, freq.nochim, track, taxa.sbering, file = "DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-SBERING.Rdata")
save(samdf, arctic.seqtab, freq.nochim, track, taxa.arctic, file = "DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-ARCTIC.Rdata")

