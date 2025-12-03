library(patchwork)

# Map: keep legend for Predator
final_map2 <- final_map +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    legend.position = "bottom",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12)
  )

# ANOVA theme for PowerPoint
anova_theme <- theme(
  text         = element_text(size = 14),
  axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
  axis.text.y  = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  plot.title   = element_text(size = 14, hjust = 0.5)
)

P_Marker   <- P_Marker   + anova_theme
P_Location <- P_Location + anova_theme
P_Sex      <- P_Sex      + anova_theme
P_Season   <- P_Season   + anova_theme
P_LocBody  <- P_LocBody  + anova_theme
P_Predator <- P_Predator + anova_theme
P_Year     <- P_Year     + anova_theme

anova_panel <- (P_Marker + P_Location + P_Sex) /
  (P_Season + P_LocBody + P_Predator) /
  P_Year

# Final layout: Map | ANOVAs
final_plot2 <- anova_panel

final_plot2

ggsave(
  "final_plot2_12s_noNMDS.png",
  plot   = final_plot2,
  width  = 12,
  height = 8,
  units  = "in",
  dpi    = 300,
  bg     = "white"
)

