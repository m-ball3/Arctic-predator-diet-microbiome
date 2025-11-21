library(dplyr)
library(stringr)
library(ggplot2)

# 1. Read your data
samdf <- read.csv("metadata/ADFG_dDNA_sample_metadata.csv")
samdf$Species <- tolower(samdf$Species)
labdf <- read.csv("metadata/ADFG_dDNA_labwork_metadata.csv")
combined_df <- read.csv("Deliverables/allrows-abundance.csv")

# 2. Create base specimen IDs to aid matching
labdf <- labdf %>%
  mutate(base_id = str_remove(Specimen.ID, "-[SF]$"))

# 3. Find specimens with both stomach and fecal samples
base_ids_with_both <- labdf %>%
  group_by(base_id) %>%
  summarize(has_stomach = any(Location.in.body == "stomach"),
            has_feces = any(Location.in.body == "feces")) %>%
  filter(has_stomach & has_feces) %>%
  pull(base_id)

# 4. Filter labdf and abundance for those IDs
labdf_filtered <- labdf %>%
  filter(base_id %in% base_ids_with_both) %>%
  select(Specimen.ID, Location.in.body, LabID)

# 5. Remove replicates by LabID
replicates_to_remove <- c("WADE-003-124", "WADE-003-115")
labdf_filtered <- labdf_filtered %>%
  filter(!LabID %in% replicates_to_remove)

# 5. Join lab metadata with abundance and taxonomy data
result <- labdf_filtered %>%
  left_join(combined_df %>% 
              select(Specimen.ID, Abundance, Predator, Location, Month, Year, Sex, 
                     Kingdom, Phylum, Order, Class, Family, Genus.y, Species.y, Marker), 
            by = "Specimen.ID")



# 7. OPTIONAL: Filter out rows without abundance for downstream quantitative analyses
result <- result %>% filter(!is.na(Abundance))

# 8. Calculate relative abundance per sample (Specimen.ID)
result <- result %>%
  group_by(Specimen.ID, Marker) %>%
  mutate(rel_abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

result <- result %>%
  mutate(Specimen.ID = str_trim(Specimen.ID),
         base_id = str_remove(Specimen.ID, "-[SF]$"))

# 9. Plot: Relative abundance by sample and location
rel.plot <- ggplot(result, aes(x = base_id, y = rel_abundance, fill = Species.y)) +
  geom_col(position = "stack") +
  facet_grid(Location.in.body ~ Marker, switch = "y") +  # Facet rows by Location.in.body, columns by marker
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "bottom")

print(rel.plot)

ggsave("Deliverables/Comparisons/stomach.vs.fecal.png", plot = rel.plot, width = 16, height = 8, units = "in", dpi = 300)


# renames results table because I'll want to use this downstream?

# Saves combined_df as a .csv for Amy's abundance table
write.csv(combined_df, "Deliverables/allrows-abundance.csv")


