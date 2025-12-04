# Final Plot 2

# MAP
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
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Alaska Sample Sites", x = NULL, y = NULL)
final_map


#NMDS
library(dplyr)

plot_data <- plot_data %>%
  mutate(
    Season = case_when(
      Month %in% c("DEC", "JAN", "FEB") ~ "Winter",
      Month %in% c("MAR", "APR", "MAY") ~ "Spring",
      Month %in% c("JUN", "JUL", "AUG") ~ "Summer",
      Month %in% c("SEP", "OCT", "NOV") ~ "Fall",
      TRUE                              ~ NA_character_
    )
  )
nmds.12.season <- ggplot(
  plot_data,
  aes(x = MDS1, y = MDS2, color = Predator, shape = Season)
) +
  stat_ellipse(
    aes(fill = Predator),
    geom  = "polygon",
    alpha = 0.25,
    color = NA,
    level = 0.95,
    type  = "t"
  ) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Bray NMDS",
    x = "NMDS1",
    y = "NMDS2",
    color = "Predator",
    fill  = "Predator",
    shape = "Season"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16)
  )

nmds.12.season


# ------------------------------------------------------------------
# Makes a large plot with one way ANOVAs (no color meaning)
# ------------------------------------------------------------------

custom_titles <- c(
  Marker            = "Genetic Marker",
  Location          = "Location in Alaska",
  Sex               = "Sex",
  Season            = "Season",
  `Location.in.body`= "Location in Body",
  Predator          = "Predator",
  Year              = "Year"
)

plot_anova_letters <- function(df, factor_col, response = "Species.Richness",
                               show_yaxis = TRUE) {
  factor_col <- enquo(factor_col)
  varname <- quo_name(factor_col)
  
  if (varname == "Location.in.body") {
    df <- df %>%
      filter(Location.in.body != "" &
               !is.na(Location.in.body) &
               Location.in.body %in% c("stomach", "feces"))
  }
  
  df[[varname]] <- as.factor(df[[varname]])
  
  # If not enough levels, just boxplot
  if (nlevels(df[[varname]]) < 2) {
    message("Not enough levels in ", varname, " for ANOVA. Plotting boxplot only.")
    
    p <- ggplot(df, aes_string(x = varname, y = response)) +
      geom_boxplot(fill = "grey80", color = "grey30") +
      theme_light() +
      theme(
        legend.position   = "none",
        axis.text.x       = element_text(angle = 45, hjust = 1),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank()
      ) +
      xlab(NULL) +
      ggtitle(varname)
    
    if (show_yaxis) {
      p <- p + ylab(response)
    } else {
      p <- p + ylab(NULL) +
        theme(axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_blank())
    }
    return(p)
  }
  
  # ANOVA + Tukey letters
  formula <- as.formula(paste(response, "~", varname))
  aov_res <- aov(formula, data = df)
  
  tukey <- tryCatch({
    agricolae::HSD.test(aov_res, varname, group = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(tukey) || is.null(tukey$groups)) {
    message("No Tukey groups for ", varname, ". Plotting without letters.")
    letters_df <- NULL
  } else {
    letters_df <- data.frame(
      Level  = rownames(tukey$groups),
      Letter = tukey$groups$groups
    )
  }
  
  means <- aggregate(as.formula(paste(response, "~", varname)),
                     data = df, FUN = mean)
  means <- means[order(means[[varname]]), ]
  if (!is.null(letters_df)) {
    means <- merge(means, letters_df,
                   by.x = varname, by.y = "Level", all.x = TRUE)
  }
  
  p <- ggplot(df, aes_string(x = varname, y = response)) +
    geom_boxplot(fill = "grey80", color = "grey30") +
    theme_light() +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(angle = 45, hjust = 1),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(hjust = 0.5)
    ) +
    xlab(NULL) +
    ggtitle(custom_titles[[varname]])
  
  if (show_yaxis) {
    p <- p + ylab(response)
  } else {
    p <- p + ylab(NULL) +
      theme(axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank())
  }
  
  if (!is.null(letters_df)) {
    p <- p +
      geom_text(
        data = means,
        aes_string(
          x     = varname,
          y     = max(df[[response]], na.rm = TRUE) * 1.05,
          label = "Letter"
        ),
        inherit.aes = FALSE,
        size = 5
      )
  }
  
  return(p)
}

# Make individual panels
P_Marker   <- plot_anova_letters(species_richness_long_all, Marker,           show_yaxis = TRUE)
P_Location <- plot_anova_letters(species_richness_long_all, Location,         show_yaxis = FALSE)
P_Sex      <- plot_anova_letters(species_richness_long_all, Sex,              show_yaxis = FALSE)
P_Season   <- plot_anova_letters(species_richness_long_all, Season,           show_yaxis = TRUE)
P_LocBody  <- plot_anova_letters(species_richness_long_all, Location.in.body, show_yaxis = FALSE)
P_Predator <- plot_anova_letters(species_richness_long_all, Predator,         show_yaxis = FALSE)
P_Year     <- plot_anova_letters(species_richness_long_all, Year,             show_yaxis = TRUE)

# Combine with patchwork
combined_plot <- (P_Marker + P_Location + P_Sex) /
  (P_Season + P_LocBody + P_Predator) /
  P_Year

combined_plot


(final_map / nmds.12.season) | combined_plot
