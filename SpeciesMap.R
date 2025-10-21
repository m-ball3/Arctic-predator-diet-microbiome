library(ggplot2)
library(tidyverse)

# Loads df
samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")
# Gets sample metadata; filters out NA's (shipment 1)
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")%>%
  filter(!is.na(LabID))

# Renames "species" column to "Predator"
samdf <- dplyr::rename(samdf, Predator = Species)

# Creates a column corresponding ADFG sample IDs with WADE sample IDs
samdf <- samdf %>%
  left_join(
    labdf %>% select(Specimen.ID, Repeat.or.New.Specimen., LabID),
    by = c("Specimen.ID", "Repeat.or.New.Specimen.")
  )

# beluha whale == Beluga whale
samdf$Predator <- tolower(samdf$Predator)

# RemoveS rows where LabID is NA 
## because some samples were damaged in shipment and weren't extracted
samdf <- samdf[!is.na(samdf$LabID), ]

# Remove duplicate Specimen.ID, keep first occurrence (ignores different LabID)
samdf_unique <- samdf[!duplicated(samdf$Specimen.ID), ]

# Creates a plot for the predator species
ggplot(samdf_unique, aes(x = Predator)) +
  geom_bar(fill = "steelblue") +
  geom_text(stat = "count", # Use the count statistic
            aes(label = after_stat(count)), # Map the computed count to the label
            vjust = -0.5, # Adjust vertical position to be above the bar
            color = "black", # Set label color
            size = 3 # Set label size
  )+
  labs(title = "Number of Samples by Predator",
       x = "Predator Species",
       y = "Number of Samples") +
  theme_minimal()
