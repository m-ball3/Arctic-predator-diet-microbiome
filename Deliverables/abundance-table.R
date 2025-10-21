# ------------------------------------------------------------------
# SPECIES BY ABUNDANCE TABLE
# ------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(openxlsx)

# Loads in data
combined_df <- read.csv("Deliverables/allrows-abundance.csv")

# Checks relativbe abundance is calculated correctly 
equal.one <- combined_df %>%
  group_by(Sample, Marker) %>%
  summarise(total_abundance = sum(Abundance))

# Creates a table with only desired columns
abundance_df <- combined_df %>%
  select(Abundance, Sample, Marker, Species)%>%
  drop_na()

abundance_by_marker <- abundance_df %>%
  group_by(Marker, Species) %>%
  summarise(species_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Marker) %>%
  mutate(total_abundance = sum(species_abundance),
         prop_abundance = species_abundance / total_abundance) %>%
  arrange(Marker, desc(prop_abundance))

# By marker
ggplot(abundance_by_marker, aes(x = reorder(Species, prop_abundance), y = prop_abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Marker, ncol = 1, scales = "free_x") +
  labs(title = "Proportional Abundance of Species by Marker",
       x = "Species",
       y = "Proportional Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Total across markers --> issue here is to not count same samples twice!

### CREATE EXCEL
# Filter and order data for each marker
sheet_12s <- abundance_by_marker %>%
  filter(Marker == "12s") %>%
  arrange(desc(prop_abundance)) %>%
  select(Marker, Species, prop_abundance)

sheet_16s <- abundance_by_marker %>%
  filter(Marker == "16s") %>%
  arrange(desc(prop_abundance)) %>%
  select(Marker, Species, prop_abundance)

# Create workbook and add sheets
wb <- createWorkbook()
addWorksheet(wb, "12s")
addWorksheet(wb, "16s")
writeData(wb, "12s", sheet_12s)
writeData(wb, "16s", sheet_16s)

getwd()
# Save the workbook as .xlsx (multi-sheet Excel file)
saveWorkbook(wb, "Deliverables/marker_species_proportions.xlsx", overwrite = TRUE)



# old code
# tries to deal with replicates 
## replicates are not a thing here --> need to address this


# Vectors of samples that have technical replicates:
## "?" means I haven't figured it out yet --> is it possible that 
## EB22PH005-S was mistyped and it actually EB23PH005-S???

# original <- c("WADE-003-080", "WADE-003-086", "WADE-003-082", "WADE-003-085", 
#               "WADE-003-075", "WADE-003-071", "WADE-003-077", "WADE-003-084", 
#               "WADE-003-081", "?", "WADE-003-115", "WADE-003-111", 
#               "WADE-003-128", "WADE-003-137", "WADE-003-143")
# # Vectors of samples that are technical replicates (in same order as original): 
# replicates <- c("WADE-003-093", "WADE-003-094", "WADE-003-095", "WADE-003-096", 
#                 "WADE-003-097", "WADE-003-098", "WADE-003-099", "WADE-003-100-C",
#                 "WADE-003-101-C", "WADE-003-122", "WADE-003-123", "WADE-003-124", 
#                 "WADE-003-145", "WADE-003-146", "WADE-003-147")
# 
# # Vectors of samples that have a cleaned and uncleaned version
# cleaned <- c("WADE-003-101-C", "WADE-003-101-C", "WADE-003-116-C", "WADE-003-117-C",
#              "WADE-003-118-C", "WADE-003-127-C")
# 
# uncleaned <- c("WADE-003-100-UC", "WADE-003-101-UC", "WADE-003-116-UC", "WADE-003-117-UC",
#                "WADE-003-118-UC","WADE-003-127-UC")

# # Creates replicate mapping for original and replicate samples
# rep_map <- setNames(original, replicates)
# 
# # Creates cleaned-uncleaned mapping (choose which represents the group)
# clean_map <- setNames(cleaned, uncleaned)
# 
# 
# # Combine all mappings into one named vector for convenience
# # CURRENTLY EXCLUDING "?"
# combined_map <- c(rep_map, clean_map)

# # FUNCTION: Maps samples in combined_map to sample group or keep original if no mapping exists
# ## (e.g. groups uncleaned and cleaned samples into the same group
# ## to later figure out how to deal with)
# map_sample_group <- function(sample) {
#   if (sample %in% names(combined_map)) {
#     combined_map[sample]
#   } else {
#     sample
#   }
# }
# 
# # Uses function on abundance_df (to group based on replication or not)
# abundance_df <- abundance_df %>%
#   rowwise() %>%
#   mutate(SampleGroup = map_sample_group(Sample)) %>%
#   ungroup()
# 
# # Renormalizes
# 
# # Sums abundances across replicates (SampleGroup + Marker + Species)
# abundance_summed <- abundance_df %>%
#   group_by(SampleGroup, Marker, Species) %>%
#   summarise(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Re-normalizea abundances per SampleGroup and Marker
# abundance_normalized <- abundance_summed %>%
#   group_by(SampleGroup, Marker) %>%
#   mutate(Relative_Abundance = Total_Abundance / sum(Total_Abundance, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Step 4: Check sums now
# check_sums <- abundance_normalized %>%
#   group_by(SampleGroup, Marker) %>%
#   summarise(total_abundance = sum(Relative_Abundance, na.rm = TRUE)) %>%
#   ungroup()
# 
# # LOST SPECIES ASSIGNMENT, WHICH IS THE POINT LOL!! FIGURE OUT HOW THAT 
# # GETS IN HERE
# 
