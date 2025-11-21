# ------------------------------------------------------------------
# SPECIES BY ABUNDANCE TABLE
# ------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(openxlsx)
library(gridExtra)
library(RColorBrewer)

# Loads in phyloseq objects
ps.12s <- readRDS("ps.12s")
ps.16s <- readRDS("ps.16s")

# Gets the samples
sam.12s <- sample_names(ps.12s)
sam.16s <- sample_names(ps.16s)

# Gets the relative abundance of each species by marker
rel.12s <- transform_sample_counts(ps.12s, function(x) x / sum(x))
rel.16s <- transform_sample_counts(ps.16s, function(x) x / sum(x))

df_12s <- psmelt(rel.12s)
df_16s <- psmelt(rel.16s)
df_12s$Marker <- "12S"
df_16s$Marker <- "16S"

# Creates the combined df
combined_df <- bind_rows(df_12s, df_16s)

# Checks relativbe abundance is calculated correctly 
equal.one <- combined_df %>%
  group_by(Sample, Marker) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  print(n = 22)

# Saves combined_df as a .csv for Amy's abundance table
write.csv(combined_df, "Deliverables/allrows-abundance.csv")

# Checks relativbe abundance is calculated correctly 
equal.one <- combined_df %>%
  group_by(Sample, Marker, Predator) %>%
  summarise(total_abundance = sum(Abundance))

# Creates a table with only desired columns
abundance_df <- combined_df %>%
  dplyr::select(Abundance, Sample, Marker, Species.y, Predator)%>%
  drop_na()

abundance_by_marker <- abundance_df %>%
  group_by(Marker, Species.y, Predator) %>%
  summarise(species_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Marker, Predator) %>%
  mutate(total_abundance = sum(species_abundance),
         prop_abundance = species_abundance / total_abundance) %>%
  arrange(Marker, desc(prop_abundance))

# By marker
ggplot(abundance_by_marker, aes(x = reorder(Species.y, prop_abundance), y = prop_abundance, fill = Species.y)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Marker, ncol = 1, scales = "free_x") +
  labs(title = "Proportional Abundance of Species by Marker",
       x = "Species",
       y = "Proportional Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# By marker and predator
ggplot(abundance_by_marker, aes(x = reorder(Species.y, prop_abundance), y = prop_abundance, fill = Species.y)) +
  geom_bar(stat = "identity") +
  facet_grid(Predator ~ Marker, scales = "free_x") +
  labs(title = "Proportional Abundance of Species by Marker and Predator",
       x = "Species",
       y = "Proportional Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Total across markers; predator

# Split abundance_by_marker by Predator and Marker
split_data <- split(abundance_by_marker, list(abundance_by_marker$Predator, abundance_by_marker$Marker), drop = TRUE)

# Create a new workbook
wb <- createWorkbook()

# Add a sheet and write data for each Predator-Marker combination
for(name in names(split_data)) {
  # Clean sheet name by replacing "." with "_" or other characters to avoid Excel limitations
  sheet_name <- gsub("\\.", "_", name)
  
  # Select columns you want to export
  data_to_write <- split_data[[name]] %>%
    dplyr::select(Marker, Species.y, prop_abundance, Predator)
  
  # Add worksheet named by predator_marker, truncating if needed for Excel sheet name length limits (max 31 chars)
  addWorksheet(wb, sheetName = substr(sheet_name, 1, 31))
  
  # Write data
  writeData(wb, sheet = substr(sheet_name, 1, 31), x = data_to_write)
}

# Save workbook
saveWorkbook(wb, "Deliverables/predator_marker_species_proportions_FIXED.xlsx", overwrite = TRUE)












# ------------------------------------------------------------------
# PLAYING WITH PLOTS
# ------------------------------------------------------------------
# # Reads in the sheets
# file <- "Deliverables/predator_marker_species_proportions.xlsx"
# wb_sheets <- getSheetNames(file)
# 
# # Reads all sheets into a list of data frames
# list_of_dfs <- lapply(wb_sheets, function(sheet) {
#   read.xlsx(file, sheet = sheet)
# })
# 
# # Extracts all unique species present in all dfs
# all_species <- unique(unlist(lapply(list_of_dfs, function(df) df$Species.y)))
# 
# # Generates distinct colors for all species
# palette_colors <- brewer.pal(min(length(all_species), 12), "Set3")  # use 12 max colors from Set3 palette as example
# if(length(all_species) > length(palette_colors)) {
#   # Expands colors if more species than palette colors
#   palette_colors <- colorRampPalette(palette_colors)(length(all_species))
# }
# species_colors <- setNames(palette_colors, all_species)
# 
# # Applies manual color scale in plot loop
# plots <- list()
# 
# for(i in seq_along(list_of_dfs)) {
#   df <- list_of_dfs[[i]]
#   
#   # Use the sheet name from wb_sheets vector for title
#   sheet_name <- wb_sheets[i]
#   
#   df_filtered <- df %>% filter(prop_abundance > 0)
#   
#   # Skip if no data after filtering
#   if(nrow(df_filtered) == 0) next
#   
#   p <- ggplot(df_filtered, aes(x = reorder(Species.y, prop_abundance), y = prop_abundance, fill = Species.y)) +
#     geom_bar(stat = "identity") +
#     scale_fill_manual(values = species_colors) +  # consistent colors across all plots
#     labs(title = sheet_name,
#          x = "",
#          y = "Proportional Abundance") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#   
#   plots[[sheet_name]] <- p
# }
# 
# # Now only non-empty plots are in `plots`
# grid.arrange(grobs = plots, ncol = 3)
# 



# ------------------------------------------------------------------
# old code
# ------------------------------------------------------------------
# ### CREATE EXCEL
# # Filter and order data for each marker
# sheet_12s <- abundance_by_marker %>%
#   filter(Marker == "12S", Predator) %>%
#   arrange(desc(prop_abundance)) %>%
#   select(Marker, Species, prop_abundance, Predator)
# 
# sheet_16s <- abundance_by_marker %>%
#   filter(Marker == "16S") %>%
#   arrange(desc(prop_abundance)) %>%
#   select(Marker, Species, prop_abundance, Predator)
# 
# # Create workbook and add sheets
# wb <- createWorkbook()
# addWorksheet(wb, "12S")
# addWorksheet(wb, "16S")
# writeData(wb, "12S", sheet_12s)
# writeData(wb, "16S", sheet_16s)
# 
# getwd()
# # Save the workbook as .xlsx (multi-sheet Excel file)
# saveWorkbook(wb, "Deliverables/marker_species_proportions.xlsx", overwrite = TRUE)




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
