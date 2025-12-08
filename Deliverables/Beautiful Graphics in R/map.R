library(ggplot2)
library(dplyr)
library(tidygeocoder)

# Read and filter for 12S samples
combined_df <- read.csv("metadata/workable.df.csv") %>%
  filter(Marker == "12S")

# Get Alaska and Russia polygons
alaska <- subset(map_data("world"), region == "USA" & long < -140 & lat > 50)
russia <- subset(map_data("world"), region == "Russia" & long < -165 & lat > 60)

# Geocode unique locations
locations <- tibble(Location = unique(combined_df$Location)) %>%
  distinct() %>%
  geocode(Location, method = 'osm', lat = latitude, long = longitude)

# Merge coordinates with sample data, add minimal jitter if overlap
combined_df_with_coords <- combined_df %>%
  left_join(locations, by = "Location") %>%
  group_by(latitude, longitude) %>%
  mutate(
    lat_plot = ifelse(n() == 1, latitude, latitude + runif(n(), -0.01, 0.01)),
    long_plot = ifelse(n() == 1, longitude, longitude + runif(n(), -0.01, 0.01))
  ) %>%
  ungroup()

# Plot map: Alaska, Russia, points colored by Predator, labels by Location
final_map <- ggplot() +
  geom_polygon(data = alaska, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  geom_polygon(data = russia, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  geom_point(data = combined_df_with_coords,
             aes(x = long_plot, y = lat_plot, color = Predator),
             size = 6, alpha = 0.9) +
  geom_text(data = combined_df_with_coords,
            aes(x = long_plot, y = lat_plot, label = Location),
            size = 4, vjust = -1, fontface = "bold") +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(xlim = c(-180, -140), ylim = c(50, 74),
              ratio = 1 / cos(mean(alaska$lat, na.rm = TRUE) * pi / 180)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "left",
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Alaska Sample Sites", x = NULL, y = NULL)
final_map

# Save as PNG
ggsave("alaska_sample_map_no_bathy.png", final_map, width = 12, height = 9, dpi = 350, bg = "white")



