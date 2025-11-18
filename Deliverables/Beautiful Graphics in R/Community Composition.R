# ------------------------------------------------------------------
# Nice Bar Plot
# intra-ordinal community composition by Specimen.ID by Predator 
# Colours represent major prey families represented, 
# shades of each colour represent genus species binomial in each family
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Sets up Environment & Loads in Data; Cleans
# ------------------------------------------------------------------
# Loads libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(purrr)

# Loads R data file
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_16SP1+2.Rdata")

# Loads processed abundance table
otu_prop <- read.csv("./Deliverables/16S/ADFG_16s_relative_speciesxsamples.csv")

# Loads metadata
meta <- read.csv("./metadata/ADFG_dDNA_sample_metadata.csv")
# Renames "species" column to "Predator" & makes all characters lowercase
meta <- dplyr::rename(meta, Predator = Species)
meta <- meta %>%
  mutate(Predator = tolower(Predator))

#Loads taxonomic annotation table mapping species
tax <- as.data.frame(merged.taxa) # should have columns: Species.y / Family / Predator
tax <- dplyr::rename(tax, Species = Species.y)

# ------------------------------------------------------------------
# Creates a long format df with needed columns
# ------------------------------------------------------------------

# Prepare long format
otu_long <- otu_prop %>%
  pivot_longer(-Specimen.ID, names_to = "Species", values_to = "Abundance", )

# Merge with metadata to get Predator
otu_long <- otu_long %>%
  left_join(meta %>% select(Specimen.ID, Predator), by = "Specimen.ID")

# Replace dots with spaces for matching
otu_long <- otu_long %>%
  mutate(Species = gsub("\\.", " ", Species))

# Merge with taxonomy to get Family
otu_long <- otu_long %>%
  left_join(tax, by = "Species")

# Aggregate proportional abundance by sample/order/subgroup/week/year
plot_data <- otu_long %>%
  group_by(Predator, Specimen.ID, Family, Species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# For each sample, calculate abundance proportions by order/subgroup
plot_data <- plot_data %>%
  group_by(Specimen.ID) %>%
  mutate(Proportion = Abundance / sum(Abundance)) %>%
  ungroup()

# Now, aggregate by week, year, order, subgroup for plotting (mean or sum, as needed)
plot_data_sum <- plot_data %>%
  group_by(Specimen.ID, Family, Predator, Species) %>%
  summarize(Proportion = mean(Proportion), .groups = "drop")

# Choses only those than have an abundance >0
plot_data_sum_filtered <- plot_data_sum %>%
  filter(Proportion > 0)



facet_grid(Predator ~ Family, scales = "free_x", space = "free")

# Calculates abundance of species per family
species_order <- plot_data_sum_filtered %>%
  group_by(Family, Species) %>%
  summarize(total_abund = sum(Proportion), .groups = "drop") %>%
  arrange(Family, desc(total_abund))

species_order <- species_order %>%
  group_by(Family) %>%
  mutate(species_rank = rank(-total_abund, ties.method = "first"))


# Merges back into main table
plot_data_sum_filtered <- plot_data_sum_filtered %>%
  left_join(species_order %>% select(Family, Species, species_rank), by = c("Family", "Species"))


# ------------------------------------------------------------------
# Creates a color gradient for each unique family
# ------------------------------------------------------------------
unique_families <- na.omit(unique(plot_data_sum_filtered$Family))
species_order_no_na <- species_order %>%
  filter(!is.na(Family))

family_gradients <- map(setNames(unique_families, unique_families), function(fam) {
  n_sp <- sum(species_order_no_na$Family == fam)
  if (!is.na(n_sp) && n_sp > 0) {
    colorRampPalette(brewer.pal(9, "YlGnBu"))(n_sp)
  } else {
    character(0)
  }
})

# Assigns each species a color by rank
species_colors <- species_order_no_na %>%
  mutate(color = map2_chr(Family, species_rank, ~ family_gradients[[.x]][.y]))

# Named vector for scale_fill_manual
fill_colors <- setNames(species_colors$color, species_colors$Species)

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
ggplot(plot_data_sum_filtered, aes(x = Specimen.ID, y = Proportion, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(Predator ~ Family, scales = "free_x", space = "free") +
  scale_fill_manual(values = fill_colors) +
  theme_minimal() +
  labs(x = "Sample", y = "Proportional abundance") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "black", fill = "white", size = 1.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.spacing = unit(0.3, "lines")
  )




