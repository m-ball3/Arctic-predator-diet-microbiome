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
library(lmtest)
library(ggeffects)
library(png)
library(grid)

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

# df that includes all sucessfully amplified ssamples in 12s and 16s 
# (not just the ones that were successful for both 12s and 16s)
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
  select(Specimen.ID, Species.Richness, Sex, Predator)%>%
  mutate(Sex = as.factor(Sex))%>%
  mutate(Predator = as.factor(Predator))


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

pred.p <- ggplot(sex.df, aes(x = Predator, y = Species.Richness)) +
  geom_boxplot(fill = "lightblue", color = "grey30") +
  theme_light() +
  theme(
    legend.position   = "none",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )
pred.p

sex.p <- ggplot(sex.df, aes(x = Sex, y = Species.Richness)) +
  geom_boxplot(fill = "skyblue3", color = "grey30") +
  theme_light() +
  theme(
    legend.position   = "none",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )
sex.p

pred.p | sex.p

# ------------------------------------------------------------------
# Fits Model
# ------------------------------------------------------------------
# Models

pred.lm <- lm(Species.Richness ~ Predator, sex.df, na.action = na.fail)
pred.sex <- lm(Species.Richness ~ Predator + Sex, sex.df, na.action = na.fail)
pred.sex.int <- lm(Species.Richness ~ Predator + Sex + Predator:Sex, sex.df, na.action = na.fail)

# log transformation
sex.df <- sex.df %>%
  mutate(Species.Richness.log = log(Species.Richness))

#model with log transformation
pred.lm.log <- lm(Species.Richness.log ~ Predator, sex.df, na.action = na.fail)
pred.sex.log <- lm(Species.Richness.log ~ Predator + Sex, sex.df, na.action = na.fail)
pred.sex.int.log <- lm(Species.Richness.log ~ Predator + Sex + Predator:Sex, sex.df, na.action = na.fail)

summary(pred.lm.log)
summary(pred.sex.log)
summary(pred.sex.int.log)

# Residuals
res <- residuals(pred.sex.int.log)
plot(res)

# QQ for normality
qqnorm(res)
qqline(res)

# Shapiro test for normality
shapiro.test(res) #p-value = 0.1891

# levene test for homogeneity of variance
leveneTest(pred.lm)

# Checking for influenial outliers
plot(pred.lm.log, which = 1)  # Residuals vs Fitted

# Cook's distance
plot(pred.lm, which = 4)

# MODEL CHECKS
dredge(pred.sex.int.log)

??lrtest()
lrtest(pred.lm.log, pred.sex.log, pred.sex.int.log)

# CHOSEN MODEL
summary(pred.sex.log)
plot(pred.sex.log)

# creates predicted values

pred <- ggpredict(pred.sex.log, terms = c("Predator", "Sex"))

# Model fitness plot
f.pred <- ggplot(pred,
       aes(x = x, y = predicted, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(0.7), width = 0.25) +
  geom_point(data = sex.df, 
             aes(x = Predator, y = Species.Richness.log, color = Sex),
             inherit.aes = FALSE,
             position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15), 
             alpha = 0.7, size = 2) +
  labs(title = "Predicted log(Species Richness) by Predator and Sex",
       x = "Predator", y = "log(Species Richness)", 
       fill = "Sex", color = "Sex") +
  scale_fill_manual(values = c("F" = "#1f78b4", "M" = "lightblue")) +
  scale_color_manual(values = c("F" = "#153e5e", "M" = "lightblue4")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
f.pred

# Add Lab ID to otu.abs
write.csv(sex.df,
          "./Deliverables/481/pred.sex.df.csv", 
          row.names = FALSE)

