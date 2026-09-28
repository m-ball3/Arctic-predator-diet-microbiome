# assigns taxonomy for CO1

library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)
library(dplyr)


# loads taxonomy ref files
taxref <- "Z:/Arctic predator diet sequences/CO1/MIDORI2_UNIQ_NUC_GB268_CO1_DADA2.fasta"
speciesref <- "Z:/Arctic predator diet sequences/CO1/CO1-addspecies-DB.fasta"

# loads output
# Loads in dada2 output
load("./DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_CO1P1_output-notaxa1.Rdata")

# assigns taxonomy and species

### Assign Taxonomy
taxam <- assignTaxonomy(seqtab.nochim, taxref, tryRC = TRUE, minBoot = 95)

# Assign Species
genus.species <- assignSpecies(seqtab.nochim, speciesref)

# Ensures all taxonomic levels agree for all rows within an species assignment
taxa <- as.data.frame(taxam) %>% # saves taxa as a dataframe
  rownames_to_column("ASV") %>%  # makes ASV rownames to a column
  filter(is.na(Species)) %>% # keeps only those in Species column that are NA
  left_join(as.data.frame(genus.species) %>% # adds genus.species df to tax_table by ASV column
              rownames_to_column("ASV"),  by = "ASV") %>%  
  unite(col = addSpecies, Genus.y,Species.y, sep = " ") %>%  # combines genus.species df columns into one and renames column addspecies
  ungroup() %>%  # ungroups by ASV
  mutate(addSpecies = case_when(addSpecies == "NA NA"~NA, TRUE~addSpecies)) %>% # changes addSpecies column to R NA if NA NA characters and to what is in addSpecies column when there is a species
  mutate(addSpecies = gsub(" NA", " spp.", addSpecies)) %>%  # changes NA to spp. in addSpecies column
  mutate(.grp = ifelse(is.na(addSpecies),  paste0("NA_grp_", row_number()),  addSpecies)) %>%  # creates new column where NAs are specified by rowname (NA_grp_5)
  group_by(.grp) %>%  # groups by column .grp
  mutate(Class = if (length(unique(Class)) > 1) NA else Class) %>% # makes sure column and .grp agree; if not -> NA
  mutate(Order = if (length(unique(Order)) > 1) NA else Order) %>% 
  mutate(Family = if (length(unique(Family)) > 1) NA else Family) %>% 
  mutate(Genus.x = if (length(unique(Genus.x)) > 1) NA else Genus.x) %>% 
  ungroup() %>%  # ungroups
  dplyr::rename("Genus" = Genus.x, "Species" = addSpecies) %>%  # renames weird column names to names that make sense
  select(-Species.x) %>%  # removes Species.x column
  bind_rows(as.data.frame(taxa) %>% # makes taxa a df and adds its columns
              rownames_to_column("ASV") %>%  # adds the ASVs as a column instead of rownames
              filter(!is.na(Species))) %>%  # adds all rows back in 
  mutate(.grp = ifelse(is.na(Species),  # changes .grp column to ???
                       paste0("NA_grp_", row_number()),  Species)) %>%  
  group_by(.grp) %>%  # groups by .grp
  fill(Order, Family, Genus, .direction = "updown") %>% # fills in NAs with the previous or next non NA value (adds correct value where disagreements were)
  ungroup() %>%  # ungroups
  select(-.grp) %>%  # removes .grp column
  column_to_rownames("ASV") %>%  # puts the columns ASV back to rownames
  as.matrix() # transforms df back to matrix

### Save dataco
save(seqtab.nochim, freq.nochim, track, taxa, taxam, genus.species, out, file = "WADE003-arcticpred_dada2_QAQC_CO1P1_output-withtaxa1.Rdata")
