# Creates workable df: 

#Loads in data
combined_df <- read.csv("Deliverables/allrows-abundance.csv")
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")
  
#Remove replicates by LabID
replicates_to_remove <- c("WADE-003-124", "WADE-003-115")
combined_df <- combined_df %>%
  filter(!LabID %in% replicates_to_remove)

# Adds Location in body
# Select relevant columns from labdf
labdf_subset <- labdf %>%
  select(Specimen.ID, Location.in.body)

# Join combined_df with labdf_subset by Specimen.ID
combined_df <- combined_df %>%
  left_join(labdf_subset, by = "Specimen.ID")

# Adds season based on month
combined_df <- combined_df %>%
  mutate(Season = case_when(
    Month %in% c("DEC", "JAN", "FEB") ~ "Winter",
    Month %in% c('MAR', "APR", 'MAY') ~ "Spring",
    Month %in% c('JUN', 'JUL', 'AUG') ~ "Summer",
    Month %in% c('SEP', 'OCT', 'NOV') ~ "Fall",
    TRUE ~ NA_character_
  )
  )



# Saves master df
write.csv(combined_df, "metadata/workable.df.csv")
