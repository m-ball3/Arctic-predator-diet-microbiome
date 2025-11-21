# Seals

# Sets up environemnt
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(RColorBrewer)

# Loads in data
combined_df <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")

# Renames "species" column to "Predator"
combined_df <- dplyr::rename(combined_df, Predator = Species)
combined_df$Predator <- tolower(combined_df$Predator)

# Assign season based on month
combined_df <- combined_df %>%
  mutate(Season = case_when(
      Month %in% c("DEC", "JAN", "FEB") ~ "Winter",
      Month %in% c('MAR', "APR", 'MAY') ~ "Spring",
      Month %in% c('JUN', 'JUL', 'AUG') ~ "Summer",
      Month %in% c('SEP', 'OCT', 'NOV') ~ "Fall",
      TRUE ~ NA_character_
    )
  )

# Filters for seal species
seal_species <- c("ringed seal", "bearded seal", "spotted seal")
seal_df <- combined_df %>% filter(Predator %in% seal_species)

# Number of samples per location for each seal species
p1 <- ggplot(seal_df, aes(x = Location, fill = Predator)) +
  geom_histogram(stat = "count", position = "dodge") +
  labs(title = "", x = "Location", y = "Sample Count") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")+
  scale_fill_brewer(palette = "Paired")  # Apply the BrBG palette


# Number of samples per year for each seal species
 p2 <- ggplot(seal_df, aes(x = Year, fill = Predator)) +
  geom_histogram(stat = "count", position = "stack") +
  scale_x_continuous(breaks = 2008:2024) +  # Set tick marks for each year from 2008 to 2024
  labs(title = "", x = "Year", y = "Sample Count") + 
   theme_minimal()+ 
   theme(legend.position = "bottom") +
   scale_fill_brewer(palette = "Paired")  # Apply the BrBG palette


# Number of samples per season for each seal species
p3<- ggplot(seal_df, aes(x = Season, fill = Predator)) +
  geom_histogram(stat = "count", position = "dodge") +
  labs(title = "", x = "Season", y = "") +
  theme_minimal()+
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Paired")  # Apply the BrBG palette



seals<- (p1 + p3) / p2
seals


ggsave("Deliverables/Seals.png", plot = seals, width = 16, height = 8, units = "in", dpi = 300)
