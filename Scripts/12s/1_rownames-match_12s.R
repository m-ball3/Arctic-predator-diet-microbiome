# ------------------------------------------------------------------
# FORMATS SAMDF TO HAVE ROWNAMES THAT MATCH THE ROWNAMES IN SEQTAB.NOCHIM
# THIS IS THE FIRST STEP AFTER DADA2
# ------------------------------------------------------------------

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(stringr)
library(lubridate)

# Loads in dada2 output
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12S_output.Rdata")
track_df<- as.data.frame(track)

mean(track_df$nonchim)
range(track_df$nonchim)

write.csv(track_df, "./Deliverables/12s/12s-track.csv", row.names = TRUE)

# Gets sample metadata
## TEMPORARILY RENAMES EB22PH005-S TO EB23PH005-S -----------------------------------------------------------
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!str_ends(LabID, "-UC")) %>%# Removes rows where -UC is attached to labID (no clean vs uncleaned samples for 12s)
  mutate(LabID = gsub("-C", "", LabID)) %>%# removes -C tag
  mutate(Specimen.ID = 
           ifelse(Specimen.ID == "EB22PH005-S", "EB23PH005-S", Specimen.ID))


samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv") %>%
  dplyr::rename(Predator = Species) %>% # renames species column to predator
  mutate(Predator = tolower(Predator))  # changes capitalization to all lowercase (fixes Beluga and beluga)


# ------------------------------------------------------------------
# FORMATS METADATASHEET FOR PHYLOSEQ OBJ
# ------------------------------------------------------------------

# Removes file extensions from OTU table names
rownames(seqtab.nochim) <- gsub("-MFU_S\\d+", "", rownames(seqtab.nochim))


# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(labdf %>% 
              dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
            by = c("Specimen.ID", "Repeat.or.New.Specimen."))

# Removes rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ]


# Sets row names to LabID
rownames(samdf) <- samdf$LabID

# Only keeps rows that appear in both metadata and seq.tab 
## AKA only samples that made it through all steps 
common_ids <- intersect(rownames(samdf), rownames(seqtab.nochim)) # should be 128 but its 123
samdff <- samdf[common_ids, ]
seqtab.nochimm <- seqtab.nochim[common_ids, ]


# Checks which samples are dropped
lost_from_samdf <- setdiff(rownames(samdf), rownames(samdff))
lost_from_seqtab <- setdiff(rownames(seqtab.nochim), rownames(seqtab.nochimm))


lost_from_samdf
lost_from_seqtab #we lose -122 because it is not in samdf --> WHY??

samdf <- samdff
seqtab.nochim <- seqtab.nochimm

# Checks for identical sample rownames in both
any(duplicated(rownames(samdf)))
any(duplicated(rownames(seqtab.nochim)))

all(rownames(samdf) %in% rownames(seqtab.nochim))
all(rownames(seqtab.nochim) %in% rownames(samdf))

# Samples in metadata but not in OTU table
setdiff(rownames(samdf), rownames(seqtab.nochim))

# Samples in OTU table but not in metadata
setdiff(rownames(seqtab.nochim), rownames(samdf))

# Creates master phyloseq object
ps.raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))


### Save data
save(samdf, seqtab.nochim, taxa, track, out, freq.nochim, ps.raw, file = "./Scripts/12s/rdata/rownames-match_12s.RData")

