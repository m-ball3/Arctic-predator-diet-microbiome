# ------------------------------------------------------------------
# FINAL PROJECT
# Mollie Ball
# QSCI 482
# ------------------------------------------------------------------

# Loads Libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(car)
library(ggplot2)
library(RColorBrewer)
library(agricolae)   # for post-hoc tests and letters
library(patchwork)   # for plot layout (if not already loaded)


# ------------------------------------------------------------------
# Loads in sample data
# ------------------------------------------------------------------

# Gets data
samdf <- data.frame(read.csv(file="metadata/workable.df.csv")) 

# Reads im proportion tables
otu.abs.16s <- read.csv("./Deliverables/16S/ADFG_16s_absolute_speciesxsamples.csv")
otu.abs.12s <- read.csv("./Deliverables/12S/ADFG_12s_absolute_speciesxsamples-trunc130-4.csv")

# ------------------------------------------------------------------
# keeps only the rows where the Specimen.ID appears in samdf and in otu.abs.12s or otu.abs 16s
# ------------------------------------------------------------------

# Extract Specimen.ID vectors from the OTU data frames, adjust column name if different
specimen_12s <- otu.abs.12s$Specimen.ID
specimen_16s <- otu.abs.16s$Specimen.ID

# Filter samdf to keep rows where Specimen.ID appears in either OTU data frame
samdf_filtered <- samdf %>%
  filter(Specimen.ID %in% specimen_12s | Specimen.ID %in% specimen_16s)

samdf_unique <- samdf_filtered %>%
  select(Sample, Specimen.ID, Predator, Location.in.body, Location, Month, Year, Season, Sex, Marker) %>%
  distinct(Specimen.ID, .keep_all = TRUE)

# ------------------------------------------------------------------
# Calculates the species richness 
# ------------------------------------------------------------------
richness.16s <- apply(otu.abs.16s[, -1], 1, function(x) sum(x > 0))
names(richness.16s) <- otu.abs.16s$Specimen.ID
samdf_unique$Species_Richness.16s <- richness.16s[samdf_unique$Specimen.ID]

richness.12s <- apply(otu.abs.12s[, -1], 1, function(x) sum(x > 0))
names(richness.12s) <- otu.abs.12s$Specimen.ID
samdf_unique$Species_Richness.12s <- richness.12s[samdf_unique$Specimen.ID]

rich.12s <- as.data.frame(richness.12s)

# checks what is being counted: 
sample_row <- which(otu.abs.12s$Specimen.ID == '2023295')
colnames(otu.abs.12s)[-1][otu.abs.12s[sample_row, -1] > 0]

colnames(otu.abs.12s)

# ------------------------------------------------------------------
# Dataframe to compare species richness across markers
# ------------------------------------------------------------------

# First, keep specimens that have both 12s and 16s richness
specimens_both <- samdf_unique %>%
  filter(!is.na(Species_Richness.16s) & !is.na(Species_Richness.12s)) %>%
  select(Specimen.ID, Species_Richness.16s, Species_Richness.12s)

# Then pivot longer to get a long-form data frame with Marker and Species.Richness columns
species_richness_long <- specimens_both %>%
  pivot_longer(cols = c(Species_Richness.12s, Species_Richness.16s),
               names_to = "Marker",
               values_to = "Species.Richness") %>%
  mutate(Marker = as.character(Marker)) %>%
  mutate(Marker = dplyr::recode(Marker,
                                "Species_Richness.12s" = "12s",
                                "Species_Richness.16s" = "16s"))


# Final columns: Specimen.ID, Marker, Species.Richness
species_richness_long <- species_richness_long %>%
  select(Specimen.ID, Marker, Species.Richness)

species_richness_long <- species_richness_long %>%
  left_join(samdf_unique %>% 
              select(Specimen.ID, Location, Sex, Season, Location.in.body, Predator, Year),
            by = "Specimen.ID")

species_richness_long_all <- samdf_unique %>%
  select(Specimen.ID, Species_Richness.12s, Species_Richness.16s) %>%
  pivot_longer(cols = c(Species_Richness.12s, Species_Richness.16s),
               names_to = "Marker",
               values_to = "Species.Richness") %>%
  filter(!is.na(Species.Richness)) %>%
  mutate(Marker = case_when(
    Marker == "Species_Richness.12s" ~ "12s",
    Marker == "Species_Richness.16s" ~ "16s",
    TRUE ~ Marker
  ))


species_richness_long_all <- species_richness_long_all %>%
  left_join(samdf_unique %>% 
              select(Specimen.ID, Location, Sex, Season, Location.in.body, Predator, Year),
            by = "Specimen.ID")



# ------------------------------------------------------------------
# Descriptive statistics
# ------------------------------------------------------------------
nrow(samdf_unique) # = 51

# mean species richness
mean(species_richness_long_all$Species.Richness, na.rm = TRUE) # = 4.194444

# mean species richness per marker
mean(samdf_unique$Species_Richness.16s, na.rm = TRUE) # = 4.08
mean(samdf_unique$Species_Richness.12s, na.rm = TRUE) # = 4.255319

# range of species richness
range(species_richness_long_all$Species.Richness, na.rm = TRUE) # = 4.194444

range(samdf_unique$Species_Richness.16s, na.rm = TRUE) # 1 - 11
range(samdf_unique$Species_Richness.12s, na.rm = TRUE) # 1 - 14

var(samdf_unique$Species_Richness.16s, na.rm = TRUE) # 5.993333
var(samdf_unique$Species_Richness.12s, na.rm = TRUE) # 8.672525

sd(samdf_unique$Species_Richness.16s, na.rm = TRUE) # 2.448129
sd(samdf_unique$Species_Richness.12s, na.rm = TRUE) # 2.944915

par(mfrow = c(2, 2))

hist(na.omit(samdf_unique$Species_Richness.16s), main = "Histogram: Richness 16s", xlab = "Species Richness 16s", col = "skyblue")
hist(na.omit(samdf_unique$Species_Richness.12s), main = "Histogram: Richness 12s", xlab = "Species Richness 12s", col = "lightgreen")

boxplot(na.omit(samdf_unique$Species_Richness.16s), main = "Boxplot: Richness 16s", ylab = "Species Richness 16s", col = "skyblue")
boxplot(na.omit(samdf_unique$Species_Richness.12s), main = "Boxplot: Richness 12s", ylab = "Species Richness 12s", col = "lightgreen")

par(mfrow = c(1, 1))


# ------------------------------------------------------------------
# Explores independent variables with many anovas
# ------------------------------------------------------------------

### MANY ANOVAS
# Example for Marker
anova_marker <- aov(Species.Richness ~ Marker, data = species_richness_long)
summary(anova_marker)

# Repeat for other factors
# e.g. Location, Sex, Season, Location.in.body, Predator, Year
anova_location <- aov(Species.Richness ~ Location, data = species_richness_long)
anova_sex <- aov(Species.Richness ~ Sex, data = species_richness_long)
anova_season <- aov(Species.Richness ~ Season, data = species_richness_long)
anova_locbody <- aov(Species.Richness ~ Location.in.body, data = species_richness_long)
anova_predator <- aov(Species.Richness ~ Predator, data = species_richness_long)
anova_year <- aov(Species.Richness ~ Year, data = species_richness_long)


# ------------------------------------------------------------------
# explores further & fits glm
# ------------------------------------------------------------------

# Plots histogram of species richness
hist(species_richness_long_all$Species.Richness)

#fits linear model
species.richness_glm <- lm(Species.Richness ~ 
                             (Location + Sex + Season + 
                                Location.in.body + Predator + 
                                Year + Marker)^2, 
                           data = species_richness_long_all)

# Output of summary from model
summary(species.richness_glm)

# Analysis of deviance table for the model
anova(species.richness_glm)

# Residuals
res <- residuals(species.richness_glm)
plot(res)

# QQ for normality
qqnorm(res)
qqline(res)

# Shapiro test for normality
shapiro.test(res) #p-value = 0.0004514

# levene test for homogeneity of variance 3 NEED TO DO THIS FOR OTHER INDEPENDENT VARIABLES
leveneTest(Species.Richness ~ Marker, data = species_richness_long)
# F = 0.3133 p = 0.5788


# Checking for influenial outliers
plot(species.richness_glm, which = 1)  # Residuals vs Fitted

# Cook's distance
plot(species.richness_glm, which = 4)  # Cook's distance plot










# OLD CODE
# # ------------------------------------------------------------------
# # Makes a large plot with one way anovas
# # ------------------------------------------------------------------
# 
# custom_titles <- c(
#   Marker = "Genetic Marker",
#   Location = "Location in Alaska",
#   Sex = "Sex",
#   Season = "Season",
#   `Location.in.body` = "Location in Body",
#   Predator = "Predator",
#   Year = "Year"
# )
# 
# 
# 
# plot_anova_letters <- function(df, factor_col, response = "Species.Richness",
#                                show_yaxis = TRUE) {
#   factor_col <- enquo(factor_col)
#   varname <- quo_name(factor_col)
#   
#   if (varname == "Location.in.body") {
#     df <- df %>%
#       filter(Location.in.body != "" & !is.na(Location.in.body) & Location.in.body %in% c("stomach", "feces"))
#   }
#   
#   df[[varname]] <- as.factor(df[[varname]])
#   
#   if (nlevels(df[[varname]]) < 2) {
#     message("Not enough levels in ", varname, " for ANOVA. Plotting boxplot only.")
#     p <- ggplot(df, aes_string(x = varname, y = response, fill = varname)) +
#       geom_boxplot() +
#       scale_fill_brewer(palette = "Paired") +
#       theme_light() +
#       theme(
#         legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()
#       ) +
#       xlab(NULL) +  # remove x-axis label
#       ggtitle(varname)  # use plot title for variable name
#     
#     if (show_yaxis) {
#       p <- p + ylab(response)
#     } else {
#       p <- p + ylab(NULL) + theme(axis.title.y = element_blank(),
#                                   axis.ticks.y = element_blank(),
#                                   axis.text.y = element_blank())
#     }
#     return(p)
#   }
#   
#   formula <- as.formula(paste(response, "~", varname))
#   aov_res <- aov(formula, data = df)
#   
#   tukey <- tryCatch({
#     agricolae::HSD.test(aov_res, varname, group = TRUE)
#   }, error = function(e) NULL)
#   
#   if (is.null(tukey) || is.null(tukey$groups)) {
#     message("No Tukey groups for ", varname, ". Plotting without letters.")
#     letters_df <- NULL
#   } else {
#     letters_df <- data.frame(Level = rownames(tukey$groups),
#                              Letter = tukey$groups$groups)
#   }
#   
#   means <- aggregate(as.formula(paste(response, "~", varname)), data = df, FUN = mean)
#   means <- means[order(means[[varname]]), ]
#   if (!is.null(letters_df)) {
#     means <- merge(means, letters_df, by.x = varname, by.y = "Level", all.x = TRUE)
#   }
#   
#   p <- ggplot(df, aes_string(x = varname, y = response, fill = varname)) +
#     geom_boxplot() +
#     #scale_fill_brewer(palette = "Accent") +
#     theme_light() +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       plot.title = element_text(hjust = 0.5)   # centers title
#     ) +
#     xlab(NULL) +
#     ggtitle(custom_titles[[varname]])
#   
#   
#   if (show_yaxis) {
#     p <- p + ylab(response)
#   } else {
#     p <- p + ylab(NULL) + theme(axis.title.y = element_blank(),
#                                 axis.ticks.y = element_blank(),
#                                 axis.text.y = element_blank())
#   }
#   
#   if (!is.null(letters_df)) {
#     p <- p + geom_text(data = means,
#                        aes_string(x = varname, y = max(df[[response]], na.rm = TRUE) * 1.05, label = "Letter"),
#                        inherit.aes = FALSE, size = 5)
#   }
#   return(p)
# }
# 
# P_Marker <- plot_anova_letters(species_richness_long_all, Marker, show_yaxis = TRUE)
# P_Location <- plot_anova_letters(species_richness_long_all, Location, show_yaxis = FALSE)
# P_Sex <- plot_anova_letters(species_richness_long_all, Sex, show_yaxis = FALSE)
# P_Season <- plot_anova_letters(species_richness_long_all, Season, show_yaxis = TRUE)
# P_LocBody <- plot_anova_letters(species_richness_long_all, Location.in.body, show_yaxis = FALSE)
# P_Predator <- plot_anova_letters(species_richness_long_all, Predator, show_yaxis = FALSE)
# P_Year <- plot_anova_letters(species_richness_long_all, Year, show_yaxis = TRUE)
# 
# (P_Marker + P_Location + P_Sex) / (P_Season + P_LocBody + P_Predator) / P_Year
# 
# combined_plot <- (P_Marker + P_Location + P_Sex) / (P_Season + P_LocBody + P_Predator) / P_Year
# 
# # Save the combined plot to a file
# # ggsave("./Deliverables/Beautiful Graphics in R/species_richness_anova_plots.png", combined_plot, width = 12, height = 10, dpi = 300)


