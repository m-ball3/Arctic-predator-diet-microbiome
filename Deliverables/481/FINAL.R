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

# ------------------------------------------------------------------
# Loads in sample data
# ------------------------------------------------------------------

# Gets data
samdf <- data.frame(read.csv(file="metadata/workable.df.csv")) 

# Reads im proportion tables
otu.abs.16s <- read.csv("C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE LAB/Arctic-predator-diet-microbiome/Deliverables/16S/ADFG_16s_absolute_speciesxsamples.csv")
otu.abs.12s <- read.csv("C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE LAB/Arctic-predator-diet-microbiome/Deliverables/12S/ADFG_12s_absolute_speciesxsamples-trunc130-4.csv")

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
# Descriptive statistics
# ------------------------------------------------------------------
nrow(samdf_unique) # = 51

mean(samdf_unique$Species_Richness.16s, na.rm = TRUE) # = 4.08
mean(samdf_unique$Species_Richness.12s, na.rm = TRUE) # = 4.255319

range(samdf_unique$Species_Richness.16s, na.rm = TRUE) # 1 11
range(samdf_unique$Species_Richness.12s, na.rm = TRUE) # 1 14

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
  mutate(Marker = recode(Marker,
                         "Species_Richness.12s" = "12s",
                         "Species_Richness.16s" = "16s"))

# Final columns: Specimen.ID, Marker, Species.Richness
species_richness_long <- species_richness_long %>%
  select(Specimen.ID, Marker, Species.Richness)

# ------------------------------------------------------------------
# Assumptions tests
# ------------------------------------------------------------------

# Model
model <- aov(Species.Richness ~ Marker, data = species_richness_long)
summary(model)

# Residuals
res <- residuals(model)

# QQ for normality
qqnorm(res)
qqline(res)

# Shapiro test for normality
shapiro.test(res) #p-value = 0.0004514

# levene test for homogeneity of variance
leveneTest(Species.Richness ~ Marker, data = species_richness_long)
# F = 0.3133 p = 0.5788


# Checking for influenial outliers
plot(model, which = 1)  # Residuals vs Fitted

# Cook's distance
plot(model, which = 4)  # Cook's distance plot


# think about whether these are independent.... I dont think they are



## maybe species richness is not the best metric here
## OR maybe I should not do by Marker, since I have so few that overlap by marker (n=22)
