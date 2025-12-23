# ------------------------------------------------------------------
# EXPLOREES METADATA
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Sets up the Environment and Loads in data
# ------------------------------------------------------------------
# Sets up environemnt
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(RColorBrewer)

# Loads in data
meta <- read_csv("metadata/ADFG_dDNA_sample_metadata.csv") |>
  rename(Predator = Species) |>
  mutate(Predator = tolower(Predator))

# ------------------------------------------------------------------
# Exploration - Cleans up and adds columns to metadata
# ------------------------------------------------------------------

# Assign season based on month
meta <- meta %>%
  mutate(Season = case_when(
      Month %in% c("DEC", "JAN", "FEB") ~ "Winter",
      Month %in% c('MAR', "APR", 'MAY') ~ "Spring",
      Month %in% c('JUN', 'JUL', 'AUG') ~ "Summer",
      Month %in% c('SEP', 'OCT', 'NOV') ~ "Fall",
      TRUE ~ NA_character_
    )
  )

# Fixes Tooth_age unk and NA
meta <- meta %>%
  mutate(Tooth_age = na_if(Tooth_age, "unk"),
Tooth_age = factor(Tooth_age)
)

# ------------------------------------------------------------------
# Exploration - Plots
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

# plots sample count by location, colored by predator
p2 <- ggplot(meta, aes(x = Location, fill = Predator)) +
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
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Paired")

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

# patchwork
sampbypred <- p1 / p2 / p3 / p4
sampbypred

# saves
ggsave("Deliverables/sampbypred.png",
       plot = sampbypred, width = 14, height = 16, units = "in", dpi = 300)

# creates separate multipanel plot 

# plots sample count by Tooth_age
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

# plots sample count by sex 
p6 <- ggplot(meta, aes(x = Tooth_age, fill = Predator)) +
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


sampbypred2 <- p5 / p6
sampbypred2

ggsave("Deliverables/sampbypred2.png",
       plot = sampbypred2, width = 14, height = 16, units = "in", dpi = 300)


# Filters for seal species
seal_species <- c("ringed seal", "bearded seal", "spotted seal")
seal_df <- meta %>% filter(Predator %in% seal_species)


ggsave("Deliverables/Seals.png", plot = seals, width = 16, height = 8, units = "in", dpi = 300)
