# ------------------------------------------------------------------
# EXPLOREES METADATA
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Sets up the Environment and Loads in data
# ------------------------------------------------------------------
# Sets up environemnt
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(patchwork)
library(RColorBrewer)
library(readr)
library(tidygeocoder)
library(readxl)

# Loads in data
meta <- read_csv("metadata/ADFG_dDNA_sample_metadata.csv") %>%
  rename(Predator = Species) %>%
  mutate(
    Predator = tolower(Predator),
    Location = str_replace(Location, "^Utqia\\?vik$", "Utqiaġvik"),
    Location = str_replace(Location, "Utqiaġvik", "Barrow")
  )


# ------------------------------------------------------------------
# Exploration - Cleans up and adds columns to metadata
# ------------------------------------------------------------------

# Assign season based on month
meta <- meta %>%
  mutate(
    Season = case_when(
      Month %in% c("Spring", "Summer", "Fall", "Winter") ~ Month,
      Month %in% c("DEC", "JAN", "FEB") ~ "Winter",
      Month %in% c("MAR", "APR", "MAY") ~ "Spring",
      Month %in% c("JUN", "JUL", "AUG") ~ "Summer",
      Month %in% c("SEP", "OCT", "NOV") ~ "Fall",
      TRUE ~ NA_character_
    )
  )


# Fixes Tooth_age unk and NA
# Adds in another column for age
## differences in diet for <1 year vs >= 1 yar per ADFG 
meta <- meta %>%
  mutate(
    Tooth_age = na_if(Tooth_age, "unk"),
    Tooth_age = as.numeric(as.character(Tooth_age)),  # Convert to numeric, coercing invalid to NA
    age_group = case_when(
      is.na(Tooth_age) ~ NA_character_,
      Tooth_age < 1 ~ "nursing",
      TRUE ~ "post weaning"
    ),
    Tooth_age = factor(Tooth_age)  # Factor after all mutations
  )

# Adds a column for the DB used iin diiet analysis
meta <- meta %>%
  mutate(
    sample_DB = case_when(
      Location == "Cook Inlet" ~ "cookinletDB",
      Location %in% c("Hooper Bay", "Scammon Bay") ~ "sberingDB",
      Predator == "beluga whale" & Location == "Nome" ~ "sberingDB",
      TRUE ~ "arcticDB"
    )
  )

# adds in list of successful microbiome PCR1's
mbLIST <- read_excel("metadata/mbLIST.xlsx")

# 3. Filter samdf to only the specimens in mbLIST
samdf_microbiome_successes <- meta %>%
  semi_join(mbLIST, by = "Specimen ID")

pred_counts <- samdf_microbiome_successes %>%
  count(Predator, name = "n_samples") %>%
  arrange(desc(n_samples))

pred_counts


# ------------------------------------------------------------------
# Exploration - Bar Plots
# ------------------------------------------------------------------

# sets base text size and expands y axis to avoid cutting off counts above bars
base_size <- 20
y_expand <- expansion(mult = c(0, 0.1))  # 10% extra space at top

# plots samlpe count by predator
p1 <- ggplot(meta, aes(x = Predator, fill = Predator)) +
  geom_bar() +
  geom_text(stat = "count",
            aes(label = ..count..),
            vjust = -0.3,
            size = 5) +
  labs(x = "", y = "Sample Count") +
  scale_y_continuous(expand = y_expand) +
  theme_test(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Paired")
p1

# plots sample count by Database used, colored by predator
p2 <- ggplot(meta, aes(x = sample_DB, fill = Predator)) +
  geom_bar(position = position_dodge(width = 0.8)) +
  geom_text(stat = "count",
            position = position_dodge(width = 0.8),
            aes(label = ..count..),
            vjust = -0.3,
            size = 4) +
  labs(x = "", y = "") +
  scale_y_continuous(expand = y_expand) +
  theme_test(base_size = base_size) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Paired")
p2

# plots sample count by season, colored by predator
p3 <- ggplot(meta, aes(x = Season, fill = Predator)) +
  geom_bar(position = position_dodge(width = 0.8)) +
  geom_text(stat = "count",
            position = position_dodge(width = 0.8),
            aes(label = ..count..),
            vjust = -0.3,
            size = 4) +
  labs(x = "", y = "") +
  scale_y_continuous(expand = y_expand) +
  theme_test(base_size = base_size) +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")
p3

# plots sample count by year, colored by predator
p4 <- ggplot(meta, aes(x = Year, fill = Predator)) +
  geom_bar(position = position_dodge(width = 0.8)) +
  geom_text(stat = "count",
            position = position_dodge(width = 0.8),
            aes(label = ..count..),
            vjust = -0.3,
            size = 3.5) +
  scale_x_continuous(breaks = 2008:2024) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "", y = "") +
  theme_test(base_size = base_size) +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")
p4

# plots sample count by location, colored by predator
p5 <- ggplot(meta, aes(x = Location, fill = Predator)) +
  geom_bar(position = position_dodge(width = 0.8)) +
  geom_text(stat = "count",
            position = position_dodge(width = 0.8),
            aes(label = ..count..),
            vjust = -0.3,
            size = 4) +
  labs(x = "", y = "") +
  scale_y_continuous(expand = y_expand) +
  theme_test(base_size = base_size) +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")
p5

# patchwork
sampbypred <- p1 / p3
sampbypred

sampbypred2 <- p5 / p2
sampbypred2

# saves
ggsave("Deliverables/Metadata Exploration/sampbypred.png",
       plot = sampbypred, width = 14, height = 16, units = "in", dpi = 300)

ggsave("Deliverables/Metadata Exploration/sampbypred-loc.png",
       plot = sampbypred2, width = 20, height = 16, units = "in", dpi = 300)

# creates separate multipanel plot 

# plots sample count by sex
p5 <- ggplot(meta, aes(x = Sex, fill = Predator)) +
  geom_bar(position = position_dodge(width = 0.8)) +
  geom_text(stat = "count",
            position = position_dodge(width = 0.8),
            aes(label = ..count..),
            vjust = -0.3,
            size = 4) +
  labs(x = "", y = "") +
  scale_y_continuous(expand = y_expand) +
  theme_test(base_size = base_size) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Paired")

# plots sample count by age_group
p6 <- ggplot(meta, aes(x = age_group, fill = Predator)) +
  geom_bar(position = position_dodge(width = 0.8)) +
  geom_text(stat = "count",
            position = position_dodge(width = 0.8),
            aes(label = ..count..),
            vjust = -0.3,
            size = 3.5) +
  scale_y_continuous(expand = y_expand) +
  labs(x = "Tooth Age", y = "") +
  theme_test(base_size = base_size) +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Paired")


sampbypred2 <- p6
sampbypred2

ggsave("Deliverables/Metadata Exploration/sampbypred2.png",
       plot = sampbypred2, width = 14, height = 16, units = "in", dpi = 300)


# ------------------------------------------------------------------
# Exploration - Maps
# ------------------------------------------------------------------

# Get Alaska and Russia polygons
alaska <- subset(map_data("world"), region == "USA" & long < -140 & lat > 50)
russia <- subset(map_data("world"), region == "Russia" & long < -165 & lat > 60)

# Geocode unique locations
locations <- tibble(Location = unique(meta$Location)) %>%
  distinct() %>%
  geocode(Location, method = 'osm', lat = latitude, long = longitude)

# Merge coordinates with sample data, add minimal jitter if overlap
meta_with_coords <- meta %>%
  left_join(locations, by = "Location") %>%
  group_by(latitude, longitude) %>%
  mutate(
    lat_plot = ifelse(n() == 1, latitude, latitude + runif(n(), -0.01, 0.01)),
    long_plot = ifelse(n() == 1, longitude, longitude + runif(n(), -0.01, 0.01))
  ) %>%
  ungroup()

# Computes sample counts by Location (assuming unique sites in meta_with_coords)
site_counts <- meta %>%
  count(Location, Predator, name = "n_samples") %>%
  left_join(locations %>% select(Location, latitude, longitude), by = "Location") %>%
  left_join(
    meta_with_coords %>% select(Location, Predator, long_plot, lat_plot) %>% distinct(),
    by = c("Location", "Predator")
  ) %>%
  group_by(Location) %>%
  slice_head(n = 1) %>%  # Pick first coord set per location (or average if preferred)
  ungroup()

# Plot map: Alaska, Russia, points colored by Predator, labels by Location
final_map <- ggplot() +
  geom_polygon(data = alaska, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  geom_polygon(data = russia, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  geom_point(data = site_counts,
             aes(x = long_plot, y = lat_plot, color = Predator, size = n_samples),
             alpha = 0.9) +
  geom_text(data = site_counts,
            aes(x = long_plot, y = lat_plot, label = paste0(Location, "\n(n=", n_samples, ")")),
            size = 4, vjust = -1.2, fontface = "bold") +
  scale_color_brewer(palette = "Set2") +
  scale_size_continuous(range = c(4, 12), name = "Samples") +
  coord_fixed(xlim = c(-180, -140), ylim = c(50, 74),
              ratio = 1 / cos(mean(alaska$lat, na.rm = TRUE) * pi / 180)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "left",
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Alaska Sample Sites", x = NULL, y = NULL)
final_map


# Save as PNG
ggsave("map_with_sizes.png", final_map, width = 12, height = 9, dpi = 350, bg = "white")






















### OLD CODE 
# Filters for seal species
seal_species <- c("ringed seal", "bearded seal", "spotted seal")
seal_df <- meta %>% filter(Predator %in% seal_species)


ggsave("Deliverables/Seals.png", plot = seals, width = 16, height = 8, units = "in", dpi = 300)
