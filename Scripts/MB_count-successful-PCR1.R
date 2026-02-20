# investigates how many samples from each predator species were successful in PCR 1 for Microbiome work

#Loads libraries
library(tidyverse)

# Gets sample metadata; filters out NA's (shipment 1)
## TEMPORARILY RENAMES EB22PH005-S TO EB23PH005-S -----------------------------------------------------------

successful <- read.csv("metadata/MB-successes.csv")

samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

samdf_success <- successful %>% 
  left_join(
    read.csv("metadata/ADFG_dDNA_sample_metadata.csv") %>% 
      select(Specimen.ID, Species, Sample_type),
    by = "Specimen.ID",
    relationship = "many-to-many"
  ) %>%
  rename(Predator = Species) %>%
  mutate(Predator = tolower(Predator)) %>%
  distinct(Specimen.ID, .keep_all = TRUE)

samdf_success %>% 
  count(c(Sample_type, Predator), name = "n", sort = TRUE, .drop = FALSE)


first <- c(
  "2023281","2023295","2023296","2014122","2017198","2017206","2017215",
  "2017217","2017221","2017222","2017228","2017230","2018278","2020230",
  "2021028","2021214","2021BDL-0708Sa","2021BDL-0723A","2021BDL-0723C",
  "2022161","2022165","2022Beluga-0713-SD","DL17HB001","DL17OME001",
  "DL18OME001","DL18OME002","DL18OME003","DL18OME004","DL20OME001",
  "DL20OME002","DL20OME003","DL21OTZ007","DL22OME001","DL22OME002",
  "DL22OME003","DL22OME004","DL22OME005","DL22OTZ004","DL22SCM001",
  "DL23OME001","DL23OME003","EB21GAM031-S","EB22GAM019-S","EB23PH005-S",
  "EB23PH028-S","EB23SH008-S","EB24GAM020-S","EB24GAM027-S","EB24PH055-S",
  "EB24PH075-S","PH22SH005-S","PH22SH015-S","PH22SH043-S","PH23SH005-S",
  "PH23SH013-S","PH23SH032-S","PH24GAM014-S","EB22PH005-S","PH22SH036-S",
  "EB24PH075-S","EB21GAM031-F","EB22GAM019-F","EB23PH028-F","EB23SH008-F",
  "EB23SH009-F","EB24PH055-F","EB24PH075-F","PH22SH004-F","PH22SH005-F",
  "PH22SH015-F","PH22SH036-F","PH22SH043-F","PH23GAM007-F","PH23SH005-F",
  "PH23SH013-F","PH23SH032-F","EB23PH028-F","PH22SH015-F","PH23SH032-F"
)

second <- c(
  "2023281","2023295","2023296","2014122","2017198","2017206","2017215",
  "2017217","2017221","2017222","2017228","2017230","2018278","2020230",
  "2021028","2021214","2021BDL-0708Sa","2021BDL-0723A","2021BDL-0723C",
  "2022161","2022165","2022Beluga-0713-SD","DL17HB001","DL17OME001",
  "DL18OME001","DL18OME002","DL18OME003","DL18OME004","DL20OME001",
  "DL20OME002","DL20OME003","DL21OTZ007","DL22OME001","DL22OME002",
  "DL22OME003","DL22OME004","DL22OME005","DL22OTZ004","DL22SCM001",
  "DL23OME001","DL23OME003","EB21GAM031-S","EB22GAM019-S","EB23PH005-S",
  "EB23PH028-S","EB23SH008-S","EB24GAM020-S","EB24GAM027-S","EB24PH055-S",
  "EB24PH075-S","PH22SH005-S","PH22SH015-S","PH22SH043-S","PH23SH005-S",
  "PH23SH013-S","PH23SH032-S","PH24GAM014-S","EB22PH005-S","PH22SH036-S",
  "EB21GAM031-F","EB22GAM019-F","EB23PH028-F","EB23SH008-F","EB23SH009-F",
  "EB24PH055-F","EB24PH075-F","PH22SH004-F","PH22SH005-F","PH22SH015-F",
  "PH22SH036-F","PH22SH043-F","PH23GAM007-F","PH23SH005-F","PH23SH013-F",
  "PH23SH032-F"
)

setdiff(first, second)
length(first)
length(second)
setequal(first, second)   # TRUE means same set of IDs
setdiff(first, second)    # should be character(0)
setdiff(second, first)    # also character(0) if sets are identical

length(unique(first))
length(unique(second))

duplicated_ids <- first[duplicated(first)]
duplicated_ids

samdf_success2 <- samdf_success %>% 
  mutate(
    Stomach.Goo    = tolower(Stomach.Goo),
    Intestinal.Goo = tolower(Intestinal.Goo),
    locationinbody = case_when(
      Stomach.Goo == "yes"    ~ "stomach",
      Intestinal.Goo == "yes" ~ "intestine",
      TRUE                    ~ NA_character_
    )
  )


samdf_success2 %>% 
  count(locationinbody, name = "n", sort = TRUE, .drop = FALSE)


