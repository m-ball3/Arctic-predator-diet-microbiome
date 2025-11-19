# Seals

# Sets up environemnt
library(dplyr)
library(ggplot2)
library(lubridate)

# Loads in data
combined_df <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator"
combined_df <- dplyr::rename(combined_df, Predator = Species)
combined_df$Predator <- tolower(combined_df$Predator)

# Create a year and month field (if not present)
combined_df <- combined_df %>%
  mutate(
    Year = substr(Specimen.ID, 6, 9),   # Adjust according to your Sample format
    Month = substr(Specimen.ID, 11, 12) # Likewise, adjust as necessary
  )

# Assign season based on month
combined_df <- combined_df %>%
  mutate(
    Month = as.numeric(Month),
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "Winter",
      Month %in% c(3, 4, 5) ~ "Spring",
      Month %in% c(6, 7, 8) ~ "Summer",
      Month %in% c(9, 10, 11) ~ "Fall",
      TRUE ~ NA_character_
    )
  )

# Filters for seal species
seal_species <- c("ringed seal", "bearded seal", "spotted seal")
seal_df <- combined_df %>% filter(Predator %in% seal_species)

# Number of samples per location for each seal species
ggplot(seal_df, aes(x = Location, fill = Predator)) +
  geom_histogram(stat = "count", position = "dodge") +
  labs(title = "Number of Samples per Location (Seals)", x = "Location", y = "Sample Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Number of samples per year for each seal species
ggplot(seal_df, aes(x = Year, fill = Predator)) +
  geom_histogram(stat = "count", position = "dodge") +
  labs(title = "Number of Samples per Year (Seals)", x = "Year", y = "Sample Count")

# Number of samples per season for each seal species
ggplot(seal_df, aes(x = Season, fill = Predator)) +
  geom_histogram(stat = "count", position = "dodge") +
  labs(title = "Number of Samples per Season (Seals)", x = "Season", y = "Sample Count")
