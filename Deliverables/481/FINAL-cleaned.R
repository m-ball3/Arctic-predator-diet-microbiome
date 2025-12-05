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
library(MuMIn)

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
# Calculates the species richness **I WANT TO CHANGE TO ALPHA DIVERSITY!!! --> shannon or simposon? something else?
# ------------------------------------------------------------------
richness.16s <- apply(otu.abs.16s[, -1], 1, function(x) sum(x > 0))
names(richness.16s) <- otu.abs.16s$Specimen.ID
samdf_unique$Species_Richness.16s <- richness.16s[samdf_unique$Specimen.ID]

richness.12s <- apply(otu.abs.12s[, -1], 1, function(x) sum(x > 0))
names(richness.12s) <- otu.abs.12s$Specimen.ID
samdf_unique$Species_Richness.12s <- richness.12s[samdf_unique$Specimen.ID]

rich.12s <- as.data.frame(richness.12s)

# # checks what is being counted: 
# sample_row <- which(otu.abs.12s$Specimen.ID == '2023295')
# colnames(otu.abs.12s)[-1][otu.abs.12s[sample_row, -1] > 0]
# 
# colnames(otu.abs.12s)

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
# Histogram of years
hist(na.omit(species_richness_long_all$Year), main = NA, xlab = "Year", ylab= "Number of Samples", col = "skyblue")

# Creates a df with just sames and whether they were male or female
sex.df <- species_richness_long_all %>%
  filter(Sex %in% c("M", "F")) %>%
  select(Specimen.ID, Species.Richness, Sex, Predator, Marker)

# samdf_unique is filtered by sample ID's that passed QAQC
nrow(sex.df) # = 61

#mean species richness
mean(sex.df$Species.Richness, na.rm = TRUE) # = 4.262295

# range of species richness
range(sex.df$Species.Richness, na.rm = TRUE) # = 1 - 14

# variance of species richness
var(sex.df$Species.Richness, na.rm = TRUE) # 8.263388

# standard deviation of species richness
sd(sex.df$Species.Richness, na.rm = TRUE) # 2.874611

# histogram of species richness
hist(na.omit(sex.df$Species.Richness), main = NA, xlab = "Species Richness", ylab= "Frequency", col = "skyblue")

hist(na.omit(sex.df$Species.Richness), main = "Histogram: Richness 12s", xlab = "Species Richness 12s", col = "lightgreen")
#
# boxplot(na.omit(samdf_unique$Species_Richness.16s), main = "Boxplot: Richness 16s", ylab = "Species Richness 16s", col = "skyblue")
# boxplot(na.omit(samdf_unique$Species_Richness.12s), main = "Boxplot: Richness 12s", ylab = "Species Richness 12s", col = "lightgreen")
#
# par(mfrow = c(1, 1))




# ------------------------------------------------------------------
# Fits Model
# ------------------------------------------------------------------
## when i add more variables, how does it change AOV (just predator vs just sex)


# Models

pred.lm <- lm(Species.Richness ~ Predator, sex.df, na.action = na.fail)
pred.sex <- lm(Species.Richness ~ Predator + Sex, sex.df, na.action = na.fail)
pred.sex.int <- lm(Species.Richness ~ Predator + Sex + Predator:Sex, sex.df, na.action = na.fail)
#loc.lm <- lm(Species.Richness ~ Predator + Sex + Predator:Sex, sex.df, na.action = na.fail)
summary(loc.lm)
dredge(loc.lm)

lrtest(pred.lm, pred.sex, pred.sex.int)

# Report AICS 


library(lmtest)

aov.2 <- aov(Species.Richness ~ Predator*Sex, sex.df)
summary(aov.2)


# Residuals
res <- residuals(loc.lm)
plot(res)

# QQ for normality
qqnorm(res)
qqline(res)

# Shapiro test for normality
shapiro.test(res) #p-value = 0.1123

# levene test for homogeneity of variance 3 NEED TO DO THIS FOR OTHER INDEPENDENT VARIABLES
leveneTest(aov.2)

# Checking for influenial outliers
plot(loc.lm, which = 1)  # Residuals vs Fitted

# Cook's distance
plot(species.richness_lm.log, which = 4)  # Cook's distance plot



