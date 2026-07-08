# ------------------------------------------------------------------
# THIS IS THE FOURTH STEP AFTER DADA2
## rownames-match.r, replicates-contaminated.r,
## Phyloseq-16SALL.r, must be run before this!
# ------------------------------------------------------------------

##THIS CODE IS HERE TO CHECK ON ANY 'WEIRD' SPECIES ASSIGNMENTS WE GOT FOR 12S ADFG DIET
## WHAT IS THE PROPORTION OF C. ACROSS ALL SAMPLES?
## HOW MANY SAMPLES IS IT IN?

## AFTER THESE CHECKS, MAKE A LIST OF MINOR PREY ITEMS THAT ARE NOT IN THE MAJOR PREY ITEMS LIST
## DROP NAS

# how many samples? avg. proportion (e.g. is 2% in 2 samples, its secondary prey)

library(tidyverse)


# loads in data
otu.prop <- read.csv("./Deliverables/ALL/SRKW_relative_speciesxsamples-MAJOR.csv")

# creates new df with OTU assignments as rows
out <- otu.prop %>%
  pivot_longer(
    cols = -X,
    names_to  = "species",
    values_to = "abund"
  ) %>%
  group_by(species) %>%
  summarise(
    abund = sum(abund, na.rm = TRUE) / n(),   # n() = number of rows in this species group
    .groups = "drop"
  )

out <- out %>%
  arrange(desc(abund))

writexl::write_xlsx(out, "MAJOR-species-assignment-checks.xlsx")

getwd()
