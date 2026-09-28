# ------------------------------------------------------------------
# SEPARATES SAMPLES GEOGRAPHICALLY FOR CORRECT TAXONOMIC ASSIGNMENT BY REGION
# THIS IS THE THIRD STEP AFTER DADA2
## 1_rownames-match_12S.R and 2_replicates_contaminated_12s.R must be run before this
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
load("./Scripts/12s/rdata/replicates-contaminated_12s.RData")

# loads in regional DBs
cookinletDB <- "DADA2/Ref-DB/12S/12S_Cook-Inlet-DB.fasta"
cookinletDB.sp <- "DADA2/Ref-DB/12S/12S_Cook-Inlet-addspecies-DB.fasta"

sberingDB <- "DADA2/Ref-DB/12S/12S_S-Bering-DB.fasta"
sberingDB.sp <- "DADA2/Ref-DB/12S/12S_S-Bering-addspecies-DB.fasta"

arcticDB <- "DADA2/Ref-DB/12S/12S_Arctic-DB.fasta"
arcticDB.sp<- "DADA2/Ref-DB/12S/12S_Arctic-addspecies-DB.fasta"

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
rowcount <- rownames(seqtab.nochim) # 122 as of 6/16/26

cook.ids   <- samdf$LabID[samdf$sample_DB == "cookinletDB"]
sbering.ids <- samdf$LabID[samdf$sample_DB == "sberingDB"]
arctic.ids <- samdf$LabID[samdf$sample_DB == "arcticDB"]

cook.seqtab   <- seqtab.nochim[rownames(seqtab.nochim) %in% cook.ids, ]
sbering.seqtab <- seqtab.nochim[rownames(seqtab.nochim) %in% sbering.ids, ]
arctic.seqtab  <- seqtab.nochim[rownames(seqtab.nochim) %in% arctic.ids, ]

rowcount.cook <- rownames(cook.seqtab) # 19
rowcount.sbering <- rownames(sbering.seqtab) # 20
rowcount.arctic <- rownames(arctic.seqtab) # 83
# total = 122 as of 6/16/26

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
# NEED TO ADD IN THE BLASTING CODE THAT AMY WROTE HERE TO GET BETTER ASSIGNMENTS
# ------------------------------------------------------------------



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

