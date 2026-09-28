
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


# ------------------------------------------------------------------
# Loads in Data
# ------------------------------------------------------------------
getwd()

# Loads dada2 output
load("./DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_ubiome_output.Rdata")
track_df<- as.data.frame(track)


mean(track_df$nonchim)
range(track_df$nonchim)

write.csv(track_df, "./Deliverables/ubiome/ubiome-track.csv", row.names = TRUE)

# Gets sample metadata
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!str_ends(LabID, "-UC")) %>%# Removes rows where -UC is attached to labID (no clean vs uncleaned samples for ubiome)
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
rownames(seqtab.nochim) <- sub("_S[0-9]+$", "", rownames(seqtab.nochim))

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(labdf %>% 
      dplyr::select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen."))

# Removes rows where LabID is NA (because shipment 1 was bad & thus not extracted)
samdf <- samdf[!is.na(samdf$LabID), ] #need to add in mocks, negs, and replicates to samdf

# Adds rows for negs and mock
## Specifies the IDs for mocks and negs
negs_mocks_ids <- c("mock1-782D",
                    "mock2-783D",
                    "neg1",
                    "neg2",
                    "neg2-2")

# Adds them to samdf
samdf <- samdf %>%
  mutate(Shipment = as.integer(Shipment)) %>%
  bind_rows(
    labdf %>%
      filter(LabID %in% negs_mocks_ids) %>%
      mutate(
        Shipment = NA_integer_  # NA in Shipment to integer class
      )
  )

# Sets row names to LabID
rownames(samdf) <- samdf$LabID
  
# Only keeps rows that appear in both metadata and seq.tab 
## AKA only samples that made it through all steps 
common_ids <- intersect(rownames(samdf), rownames(seqtab.nochim))

samdf <- samdf[common_ids, ]
seqtab.nochim.common <- (seqtab.nochim[common_ids, ])

lost <- setdiff(rownames(seqtab.nochim), rownames(seqtab.nochim.common))
  
lostdf <- as.data.frame(lost) # lost Arial's samples and undetermined

# after lostdf looks good, do for the real seqtab.nochim
seqtab.nochim <- (seqtab.nochim[common_ids, ])

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
                   tax_table(taxam))


### Save data
save(samdf, seqtab.nochim, taxam, track, out, freq.nochim, ps.raw, file = "./Scripts/ubiome/rdata/rownames-match_ubiome.RData")

