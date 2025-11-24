library(ggplot2)
library(dplyr)

# If you have multiple years, aggregate (summarize) by YEAR, NAME, ESR_GUILD to get maximum biomass
biomass_by_year <- biomass_stratum %>%
  group_by(YEAR, NAME, ESR_GUILD) %>%
  summarize(MAX_BIOMASS = max(MAX_BIOMASS, na.rm = TRUE))

# Create a faceted bar plot
ggplot(biomass_by_year, aes(x = factor(YEAR), y = MAX_BIOMASS, fill = ESR_GUILD)) +
  geom_col(width = 0.7, position = position_dodge()) +
  facet_wrap(~ NAME, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    panel.spacing = unit(1, "lines"),
    legend.position = "top"
  ) +
  labs(
    x = "Year",
    y = "Max Biomass",
    fill = "Guild",
    title = "Maximum Biomass by Species and Guild"
  )

ggplot(biomass_stratum, aes(x = YEAR, y = MAX_BIOMASS, color = Subregion)) +
  geom_line(size = 1.2) +
  facet_wrap(~ NAME, scales = "free_y") +
  theme_minimal() +
  labs(x = "Year", y = "Max Biomass", color = "Subregion")

ggplot(biomass_stratum, aes(x = factor(YEAR), y = MAX_BIOMASS, fill = ESR_GUILD)) +
  geom_violin() +
  facet_wrap(~ NAME, scales = "free_y") +
  theme_minimal() +
  labs(x = "Year", y = "Max Biomass", fill = "Guild")


bm_matrix <- dcast(
  biomass_stratum,
  Subregion + NAME ~ YEAR,
  value.var = "STRATUM_BIOMASS",
  fun.aggregate = max,
  na.rm = TRUE
)

bm_long <- melt(bm_matrix)
bm_long$log_value <- log10(bm_long$value + 1)  # Add 1 to avoid log(0)

p <-ggplot(bm_long, aes(x = variable, y = Subregion, fill = log_value)) +
  geom_tile() +
  facet_wrap(~ NAME) +
  scale_fill_viridis_c(option = "B", name = "Log10(Stratum Biomass + 1)") +
  theme_minimal() +
  theme(legend.position="bottom")+
  labs(x = "Year", y = "Subregion")

ggsave("stratum_biomass_heatmap.png", plot = p, width = 20, height = 7, dpi = 300)

ggplot(biomass_stratum, aes(x = YEAR, y = STRATUM_BIOMASS, color = ESR_GUILD)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ NAME) +
  theme_minimal() +
  labs(x = "Stratum Biomass", y = "Max Biomass")

