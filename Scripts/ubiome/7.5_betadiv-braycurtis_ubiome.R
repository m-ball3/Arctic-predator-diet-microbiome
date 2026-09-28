# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS ALPHA DIVERSITY METRICS ON UBIOME DATA
# THIS IS THE FIRST STATISTICAL INVESTIGATION 
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
## ALL P VALUES REPORTED ARE THE ADJUSTED P!!
# ------------------------------------------------------------------------------

# install.packages(
#   "rbiom",
#   type = "binary",
#   repos = "https://cran.rstudio.com"
# )

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(tibble)
library(stringr)
library(phyloseq)
library(car)
library(vegan)
library(rstatix)
library(Hmisc)
library(rbiom)
library(mia) # via bioconductor
library(patchwork)

getwd()
load("./Scripts/ubiome/rdata/ubiome.ps-major.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.major), "matrix")

# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------

# formatting for metadata
metadata <- data.frame(sample_data(ps.major))

class(metadata)

# optionally removes fin whale, since there is only one sample (two sub-samples of the same sample)
metadata <- metadata %>%
  dplyr::filter(Predator != "fin whale")

asv_mat <- as.data.frame(asv_mat) %>% 
  rownames_to_column(var = "row_id") %>% 
  dplyr::filter(!row_id %in% c("WADE-003-150-A", "WADE-003-150-B")) %>% 
  column_to_rownames(var = "row_id")

# formatting metadata for age class (only bearded and ringed seals)
metadata_age <- metadata %>%
  dplyr::filter(Predator %in% c("ringed seal", "bearded seal")) %>%
  dplyr::mutate(
    Age_Group = dplyr::case_when(
      Age_Group %in% c("", "Pending") ~ "Unknown",
      TRUE ~ Age_Group
    )
  )

# formatting asv_mat for age class (only bearded and ringed seals)
asv_mat_age <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

# ------------------------------------------------------------------------------
# Gets Bray-Curtis distance
# ------------------------------------------------------------------------------
# Bray–Curtis dissimilarity
bray_dist <- vegan::vegdist(
  asv_mat,
  method = "bray"
)

# Bray–Curtis dissimilarity among the seal-only age-analysis subset
bray_dist_age <- vegan::vegdist(
  asv_mat_age,
  method = "bray"
)

# Verify that metadata and distance matrix refer to the same samples
metadata <- metadata %>%
  dplyr::filter(LabID %in% labels(bray_dist)) %>%
  dplyr::arrange(match(LabID, labels(bray_dist)))

metadata_age <- metadata_age %>%
  dplyr::filter(LabID %in% labels(bray_dist_age)) %>%
  dplyr::arrange(match(LabID, labels(bray_dist_age)))


# BY PREDATOR----------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_host <- kruskal.test(Observed ~ Predator, data = metadata)
obs_test_host #  p-value = 5.066e-05

# Dunn test to determine which groups are significantly different
?dunn_test()
obs_dunn_results_host <- metadata %>%
  dunn_test(
    Observed ~ Predator,
    p.adjust.method = "hochberg"
  )
obs_dunn_results_host

# Add x- and y-coordinates for the comparison brackets
obs_dunn_plot_host <- obs_dunn_results_host %>%
  add_xy_position(x = "Predator")

# Inspect the table used for plotting
obs_dunn_plot_host

# Box-and-whisker plot, with individual samples overlaid
# BY PREDATOR
obs_plot_host <- ggplot(
  metadata,
  aes(x = Predator, y = Observed, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    obs_dunn_plot_host,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Host",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by host"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

obs_plot_host

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY PREDATOR
## bearded seal and beluga whale (adjusted p-value = 0.000827)
## ringed seal and beluga whale (adjusted p-value = 0.000571)

# gets both global and by-host summary data
summary(metadata$Observed)

metadata %>%
  group_by(Predator) %>%
  summarise(
    n = sum(!is.na(Observed)),
    median_n_Observed = median(Observed, na.rm = TRUE),
    mean_n_Observed = mean(Observed, na.rm = TRUE),
    IQR_n_Observed = IQR(Observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_n_Observed))

### biological interpretation and conclusion
# Bearded-seal and ringed-seal samples had higher observed ASV richness 
# than beluga-whale samples. There was no evidence of a difference in observed richness 
# between the two seal species



# BY SEASON---------------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_season <- kruskal.test(Observed ~ season, data = metadata)
obs_test_season # p-value = 0.8545 not significant 

# Dunn test to determine which groups are significantly different
obs_dunn_results_season <- metadata %>%
  dunn_test(
    Observed ~ season,
    p.adjust.method = "hochberg"
  )
obs_dunn_results_season

# Add x- and y-coordinates for the comparison brackets
obs_dunn_plot_season <- obs_dunn_results_season %>%
  add_xy_position(x = "season")

# Inspect the table used for plotting
obs_dunn_plot_season

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY SEASON
## no significantly differences in observed richness by season

# Box-and-whisker plot, with individual samples overlaid
# BY SEASON
obs_plot_season <- ggplot(
  metadata,
  aes(x = season, y = Observed, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    obs_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "right")

obs_plot_season

# gets both global and by-season summary data
summary(metadata$Observed)

metadata %>%
  group_by(season) %>%
  summarise(
    n = sum(!is.na(Observed)),
    median_n_Observed = median(Observed, na.rm = TRUE),
    mean_n_Observed = mean(Observed, na.rm = TRUE),
    IQR_n_Observed = IQR(Observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_n_Observed))

### biological interpretation and conclusion
#There was no evidence of a difference in observed richness 
# across seasons

## there may be some evidence (based on the plot) 
## for seasonal differences per host though

# BY LOCATION------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_location <- kruskal.test(Observed ~ Location, data = metadata)
obs_test_location # p-value = 0.001953 significant result
# at least one location is significantly different in observed richness

# Dunn test to determine which groups are significantly different
obs_dunn_results_location <- metadata %>%
  dunn_test(
    Observed ~ Location,
    p.adjust.method = "hochberg"
  )
obs_location <- as.data.frame(obs_dunn_results_location)

obs_dunn_plot_location <- obs_dunn_results_location %>%
  add_xy_position(x = "Location")

# Inspect the table used for plotting
obs_dunn_plot_location

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY LOCATION
# Cook Inlet is significantly different than Gambell (p = 0.02093657)

# Box-and-whisker plot, with individual samples overlaid
# BY LOCATION
obs_plot_location <- ggplot(
  metadata,
  aes(x = Location, y = Observed, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  )+
  stat_pvalue_manual(
    obs_dunn_plot_location,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Location",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by location"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

obs_plot_location

# gets both global and by-location summary data
summary(metadata$Observed)

metadata %>%
  group_by(Location) %>%
  summarise(
    n = sum(!is.na(Observed)),
    median_n_Observed = median(Observed, na.rm = TRUE),
    mean_n_Observed = mean(Observed, na.rm = TRUE),
    IQR_n_Observed = IQR(Observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_n_Observed))

### biological interpretation and conclusion
# Cook Inlet was significantly different than Gambel, 
# but this was likely driven by differences in host: 
# only belugas in Cook Inlet and only ringed and bearded seals in Gambel
## we already know that host is a driver of differences 

# BY LOCALE-----------------------------------------------------------------------
# Kruskal-Wallis test
obs_test_locale <- kruskal.test(Observed ~ Locale, data = metadata)
obs_test_locale # p-value = 0.0008742 significant

# Dunn test to determine which groups are significantly different
obs_dunn_results_locale <- metadata %>%
  dunn_test(
    Observed ~ Locale,
    p.adjust.method = "hochberg"
  )
obs_dunn_results_locale

obs_dunn_plot_locale <- obs_dunn_results_locale %>%
  add_xy_position(x = "Locale")

# Inspect the table used for plotting
obs_dunn_plot_locale

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY LOCALE
# Arctic is significantly different from Cook Inlet (p = 0.00222)
# Arctic is significantly different from South Bering ( p = 0.0482)

# Box-and-whisker plot, with individual samples overlaid
# BY LOCALE
obs_plot_locale <- ggplot(
  metadata,
  aes(x = Locale, y = Observed, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    obs_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Locale",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

obs_plot_locale


# gets both global and by-locale summary data
summary(metadata$Observed)

metadata %>%
  group_by(Locale) %>%
  summarise(
    n = sum(!is.na(Observed)),
    median_n_Observed = median(Observed, na.rm = TRUE),
    mean_n_Observed = mean(Observed, na.rm = TRUE),
    IQR_n_Observed = IQR(Observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_n_Observed))

### biological interpretation and conclusion
# Arctic is significantly different from Cook Inlet and South Bering
# but this was likely driven by differences in host: 
# only belugas in Cook Inlet and no ringed seals in South Bering
## we already know that host is a driver of differences 


# BY AGE--------------------------------------------------------------------------
# Kruskal-Wallis test
obs_test_age <- kruskal.test(Observed ~ Age_Group, data = metadata_age)
obs_test_age # p-value = 0.2555 not significant

# Dunn test to determine which groups are significantly different
obs_dunn_results_age <- metadata_age %>%
  dunn_test(
    Observed ~ Age_Group,
    p.adjust.method = "hochberg"
  )
obs_dunn_results_age

obs_dunn_plot_age <- obs_dunn_results_age %>%
  add_xy_position(x = "Age_Group")

# Inspect the table used for plotting
obs_dunn_plot_age

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY AGE
# none

# Box-and-whisker plot, with individual samples overlaid
# BY AGE GROUP
obs_plot_age <- ggplot(
  metadata_age,
  aes(x = Age_Group, y = Observed, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    obs_dunn_plot_age,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Age_Group",
    y = "Number of ASVs per sample",
    title = "Observed ASV richness by Age Group"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

obs_plot_age


# gets both global and by-age summary data
summary(metadata$Observed)

metadata %>%
  group_by(Age_Group) %>%
  summarise(
    n = sum(!is.na(Observed)),
    median_n_Observed = median(Observed, na.rm = TRUE),
    mean_n_Observed = mean(Observed, na.rm = TRUE),
    IQR_n_Observed = IQR(Observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_n_Observed))

### biological interpretation and conclusion
# There are no significant differences in observed richness across age groups
# in ringed and bearded seals.

# RICHNESS PLOTS ---------------------------------------------------------------

obs_plots <- (obs_plot_host + obs_plot_age + obs_plot_locale) / obs_plot_season / obs_plot_location  &
  theme(
    text = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

ggsave("Deliverables/ubiome/alphadiv/OBSRICH_ubiome-majorclass.png", plot = obs_plots, width = 40, height = 40, units = "in", dpi = 300)



# ------------------------------------------------------------------------------
# Dominance/Evenness
# ------------------------------------------------------------------------------

# Berger-Parker is a ratio between the abundance of the most abundant ASV and the 
# total abundance. 
# It is the within-sample relative abundance of the single most abundant ASV.
# Higher values indicate stronger dominance and lower community evenness.

# Compared with other dominance formula metrics, it is easy to calculate and has a 
# clear biological interpretation.

# test for normalcy
bp.norm <- shapiro.test(metadata$Berger_Parker)
bp.norm # p-value = 2.302e-05

# test homogeneity of variance
?leveneTest()
bp.var <- leveneTest(Berger_Parker ~ Predator, metadata)
bp.var # p-value = 0.924

# Does not meet assumptions for parametric tests

# tests if read depth is correlated with Berger-Parker ratios
cor.test(
  metadata$Berger_Parker,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
) 

# p-value = 0.7242; though there is a slight negative trend, there is no evidence
# that the read depth is highly correlated with the Berger-Parker ratio
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)

# BY HOST---------------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_host <- kruskal.test(Berger_Parker ~ Predator, data = metadata)
bp_test_host #  p-value = 0.0002138 significant 

# Dunn test to determine which groups are significantly different
bp_dunn_results_host <- metadata %>%
  dunn_test(
    Berger_Parker ~ Predator,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_host

### SIGNIFICANTLY DIFFERENT DOMINANCE
## bearded seal and beluga whale (adjusted p-value = 0.000366)
## bearded seal and ringed seal (adjusted p-value = 0.00206)

# Add x- and y-coordinates for the comparison brackets
bp_dunn_plot_host <- bp_dunn_results_host %>%
  add_xy_position(x = "Predator")

# Inspect the table used for plotting
bp_dunn_plot_host

# Box-and-whisker plot BY HOST, with individual samples overlaid
bp_plot_host <- ggplot(
  metadata,
  aes(x = Predator, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8)+
  stat_pvalue_manual(
    bp_dunn_plot_host,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Host",
    y = "Dominance of ASVs per sample",
    title = "ASV Dominance by host"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

bp_plot_host

# gets both global and by-host summary data by host
summary(metadata$Berger_Parker)

metadata %>%
  group_by(Predator) %>%
  summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Berger_Parker))

### biological interpretation and conclusion
# Bearded-seal samples were more frequently dominated by a single ASV than either 
# beluga-whale or ringed-seal samples



# BY SEASON----------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_season <- kruskal.test(Berger_Parker ~ season, data = metadata)
bp_test_season #  p-value = 0.05075 not significant

# Dunn test to determine which groups are significantly different
bp_dunn_results_season <- metadata %>%
  dunn_test(
    Berger_Parker ~ season,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_season

### SIGNIFICANTLY DIFFERENT DOMINANCE
## NONE

# Add x- and y-coordinates for the comparison brackets
bp_dunn_plot_season <- bp_dunn_results_season %>%
  add_xy_position(x = "season")

# Inspect the table used for plotting
bp_dunn_plot_season

# Box-and-whisker plot BY SEASON, with individual samples overlaid
bp_plot_season <- ggplot(
  metadata,
  aes(x = season, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    bp_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Host",
    y = "Dominance of ASVs per sample",
    title = "ASV Dominance by host"
  ) +
  stat_pvalue_manual(
    bp_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Dominance of ASVs per sample",
    title = "ASV Dominance by Season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "right")

bp_plot_season

# gets both global and by-host summary data by season
summary(metadata$Berger_Parker)

metadata %>%
  group_by(season) %>%
  summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Berger_Parker))

### biological interpretation and conclusion
# there are no significant differences in dominance across seasons, 
# though I bet that there are host differences across seasons (esp beluga whales)
## though I wonder if that is driven by decomp differences in hot weather vs cooler weather???

# BY LOCATION-----------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_location <- kruskal.test(Berger_Parker ~ Location, data = metadata)
bp_test_location #  p-value = 0.01363 significant
# at least one location is significantly different in dominance than the others

# Dunn test to determine which groups are significantly different
bp_dunn_results_location <- metadata %>%
  dunn_test(
    Berger_Parker ~ Location,
    p.adjust.method = "hochberg"
  )
bp_location <- as.data.frame(bp_dunn_results_location)


### SIGNIFICANTLY DIFFERENT DOMINANCE
## Point Hope is significantly different from Nome (p = 0.002701786)

# Add x- and y-coordinates for the comparison brackets
bp_dunn_plot_location <- bp_dunn_results_location %>%
  add_xy_position(x = "Location")

# Inspect the table used for plotting
bp_dunn_plot_location

# Box-and-whisker plot BY LOCATION, with individual samples overlaid
bp_plot_location <- ggplot(
  metadata,
  aes(x = Location, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    bp_dunn_plot_location,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Location",
    y = "Dominance of ASVs per sample",
    title = "ASV Dominance by location"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

bp_plot_location

# gets both global and by-host summary data by season
summary(metadata$Berger_Parker)

metadata %>%
  group_by(Location) %>%
  summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Berger_Parker))

### biological interpretation and conclusion
# Point Hope and Nome have significantly different domincane
# but this is likely driven by differences in Host
## Nome is all belugas, and Point Hope is all ringed and bearded seals

# BY LOCALE---------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_locale <- kruskal.test(Berger_Parker ~ Locale, data = metadata)
bp_test_locale #  p-value = 0.0002373 significant 
# at least one locale is significantly different from the others

# Dunn test to determine which groups are significantly different
bp_dunn_results_locale <- metadata %>%
  dunn_test(
    Berger_Parker ~ Locale,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_locale

### SIGNIFICANTLY DIFFERENT DOMINANCE
## Arctic and South Bering (p = 0.000134)
## Cook Inlet and South Bering (p = 0.0310)

# Add x- and y-coordinates for the comparison brackets
bp_dunn_plot_locale <- bp_dunn_results_locale %>%
  add_xy_position(x = "Locale")

# Inspect the table used for plotting
bp_dunn_plot_locale

# Box-and-whisker plot BY LOCALE, with individual samples overlaid
bp_plot_locale <- ggplot(
  metadata,
  aes(x = Locale, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    bp_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Locale",
    y = "Dominance of ASVs per sample",
    title = "ASV Dominance by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

bp_plot_locale

# gets both global and by-host summary data by season
summary(metadata$Berger_Parker)

metadata %>%
  group_by(Locale) %>%
  summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Berger_Parker))

### biological interpretation and conclusion
# Arctic and Cook Inlet are significantly different than South Bering, 
# but, again, this may be driven in part by Host (esp belugas)

# BY AGE GROUP----------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_age <- kruskal.test(Berger_Parker ~ Age_Group, data = metadata_age)
bp_test_age #  p-value = 0.8138 not significant

# Dunn test to determine which groups are significantly different
bp_dunn_results_age <- metadata_age %>%
  dunn_test(
    Berger_Parker ~ Age_Group,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_age

### SIGNIFICANTLY DIFFERENT DOMINANCE
## There are no significant differences by age group for ringed and bearded seals

# Add x- and y-coordinates for the comparison brackets
bp_dunn_plot_age <- bp_dunn_results_age %>%
  add_xy_position(x = "Age_Group")

# Inspect the table used for plotting
bp_dunn_plot_age

# Box-and-whisker plot BY AGE, with individual samples overlaid
bp_plot_age <- ggplot(
  metadata_age,
  aes(x = Age_Group, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    bp_dunn_plot_age,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Age Group",
    y = "Dominance of ASVs per sample",
    title = "ASV Dominance by age group"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

bp_plot_age

# gets both global and by-host summary data by season
summary(metadata_age$Berger_Parker)

metadata_age %>%
  group_by(Age_Group) %>%
  summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Berger_Parker))

### biological interpretation and conclusion
# There are taxa that are more dominant in pups vs non-pups in ringed or bearded seals
## though this may change if Host is accountded for!!!
# it looks like ringed seals may have a bit of a difference 

# DOMINANCE PLOTS --------------------------------------------------------------

bp_plots <- (bp_plot_host + bp_plot_age + bp_plot_locale) / bp_plot_season / bp_plot_location  &
  theme(
    text = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

ggsave("Deliverables/ubiome/alphadiv/BPDOM_ubiome-majorclass.png", plot = bp_plots, width = 40, height = 40, units = "in", dpi = 300)


# DOMINANT ASV -----------------------------------------------------------------

# gets the dominant ASV and the genus that it relates to
tax_df <- as.data.frame(tax_table(ps.major)) %>%
  rownames_to_column("ASV")

dominant_idx <- max.col(asv_mat, ties.method = "first")

dominant_asv_per_sample <- tibble(
  LabID = rownames(asv_mat),
  ASV = colnames(asv_mat)[dominant_idx],
  n_reads_asv = asv_mat[cbind(seq_len(nrow(asv_mat)), dominant_idx)],
  sample_reads = rowSums(asv_mat)
) %>%
  dplyr::mutate(
    dominance = n_reads_asv / sample_reads
  ) %>%
  dplyr::left_join(tax_df, by = "ASV") %>%
  dplyr::left_join(
    metadata %>% dplyr::select(LabID, Predator, season, Location, Locale, Age_Group),
    by = "LabID"
  ) %>%
  dplyr::select(
    LabID, Predator, season, Location, Locale, Age_Group, ASV,
    n_reads_asv, sample_reads, dominance,
    Kingdom, Phylum, Class
  ) %>%
  arrange(desc(dominance))

dominant_asv_per_sample

## WOULD LOVE TO MAKE A PLOT FOR DOMINANT ASV PER SAMPLE

# makes a nice summary table
dominant_asv_summary <- dominant_asv_per_sample %>%
  group_by(ASV, Kingdom, Phylum, Class) %>%
  summarise(
    n_samples_dominant = n(),
    mean_dominance = mean(dominance),
    median_dominance = median(dominance),
    max_dominance = max(dominance),
    total_reads_when_dominant = sum(n_reads_asv),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples_dominant), desc(mean_dominance))

dominant_asv_summary




# ------------------------------------------------------------------------------
# Phylogenetics
# ------------------------------------------------------------------------------

# Faith

### biological interpretation and conclusion

# ------------------------------------------------------------------------------
# Information
# ------------------------------------------------------------------------------

# Shannon

# test for normalcy
shan.norm <- shapiro.test(metadata$Shannon)
shan.norm # p-value = 3.253e-06

# test homogeneity of variance
shan.var <- leveneTest(Shannon ~ Predator, metadata)
shan.var # p-value = 0.913

# Does not meet assumptions for parametric tests

# tests if read depth is correlated with Berger-Parker ratios
cor.test(
  metadata$Shannon,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
) 

# p-value = 0.9659; though there is a slight positive trend, there is no evidence
# that the read depth is highly correlated with the Shannon index
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)

# Kruskall-Wallis test ---------------------------------------------------------

# BY HOST-------------------------------------------------------------------------------
shan_test_host <- kruskal.test(Shannon ~ Predator, data = metadata)
shan_test_host #  p-value = 0.00944 significant
# at least one host is significantly different than the others

# Dunn test to determine which groups are significantly different
shan_dunn_results_host <- metadata %>%
  dunn_test(
    Shannon ~ Predator,
    p.adjust.method = "hochberg"
  )
shan_dunn_results_host

### SIGNIFICANTLY DIFFERENT SHANNON ENTROPY
## bearded seals are significantly different from beluga whales
## bearded seals are significantly different from ringed seals
## ringed seals are not significantly different from beluga whales

# Add x- and y-coordinates for the comparison brackets
shan_dunn_plot_host <- shan_dunn_results_host %>%
  add_xy_position(x = "Predator")

# Inspect the table used for plotting
shan_dunn_plot_host

# Box-and-whisker plot, with individual samples overlaid
shan_plot_host <- ggplot(
  metadata,
  aes(x = Predator, y = Shannon, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    shan_dunn_plot_host,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Host",
    y = "Shannon Index per sample",
    title = "Shannon Index by host"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

shan_plot_host

# gets both global and by-host summary data
summary(metadata$Shannon)

metadata %>%
  group_by(Predator) %>%
  summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_shannon))

### biological interpretation and conclusion
# Shannon diversity differed across host groups in the omnibus test, 
# and three is individual host-pair contrast after Holm correction. 
# There iscorrected pairwise evidence for a 
# specific difference in Shannon diversity between bearded seals and beluga whales and ringed seals

# interpretation of all metrics collectively

# BY SEASON----------------------------------------------------------------------------
shan_test_season <- kruskal.test(Shannon ~ season, data = metadata)
shan_test_season #  p-value = 0.1109 not significant 

# Dunn test to determine which groups are significantly different
shan_dunn_results_season <- metadata %>%
  dunn_test(
    Shannon ~ season,
    p.adjust.method = "hochberg"
  )
shan_dunn_results_season

### SIGNIFICANTLY DIFFERENT SHANNON ENTROPY
## none are significant 
## there is no signficant difference in shannon entropy across seasons

# Add x- and y-coordinates for the comparison brackets
shan_dunn_plot_season <- shan_dunn_results_season %>%
  add_xy_position(x = "Season")

# Inspect the table used for plotting
shan_dunn_plot_season

# Box-and-whisker plot, with individual samples overlaid
shan_plot_season <- ggplot(
  metadata,
  aes(x = season, y = Shannon, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    shan_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Shannon Index per sample",
    title = "Shannon Index by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "right")

shan_plot_season

# gets both global and by-host summary data
summary(metadata$Shannon)

metadata %>%
  group_by(season) %>%
  summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_shannon))

### biological interpretation and conclusion
# there are no significant differences, though this may change if season by host was looked at
## there seems to be visual evidence for significant differences between beluga whale 
## microbiomes in autumn vs summer vs winter, for example. Possibly even bearded seals. 

# BY LOCATION---------------------------------------------------------------------------------
shan_test_location <- kruskal.test(Shannon ~ Location, data = metadata)
shan_test_location #  p-value = 0.1913 not significant 

# Dunn test to determine which groups are significantly different
shan_dunn_results_location <- metadata %>%
  dunn_test(
    Shannon ~ Location,
    p.adjust.method = "hochberg"
  )
shan_dunn_results_location

### SIGNIFICANTLY DIFFERENT SHANNON ENTROPY
## none are significant 

# Add x- and y-coordinates for the comparison brackets
shan_dunn_plot_location <- shan_dunn_results_location %>%
  add_xy_position(x = "Location")

# Inspect the table used for plotting
shan_dunn_plot_location

# Box-and-whisker plot, with individual samples overlaid
shan_plot_location <- ggplot(
  metadata,
  aes(x = Location, y = Shannon, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  )+
  stat_pvalue_manual(
    shan_dunn_plot_location,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Location",
    y = "Shannon Index per sample",
    title = "Shannon Index by location"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

shan_plot_location

# gets both global and by-host summary data
summary(metadata$Shannon)

metadata %>%
  group_by(Location) %>%
  summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_shannon))

### biological interpretation and conclusion
# There are no significant differences in shannon entropy across locations
## though this would probably change if location by Host was examined. 
## visually, it looks like there could be significant differences across locations for 
## beluga whales at the least. Maybe bearded seals. 

# BY LOCALE---------------------------------------------------------------------------
shan_test_locale <- kruskal.test(Shannon ~ Locale, data = metadata)
shan_test_locale #  p-value = 0.009257 significant
# at least one locale is significantly different than the others

# Dunn test to determine which groups are significantly different
shan_dunn_results_locale <- metadata %>%
  dunn_test(
    Shannon ~ Locale,
    p.adjust.method = "hochberg"
  )
shan_dunn_results_locale

### SIGNIFICANTLY DIFFERENT SHANNON ENTROPY
## Arctic is significantly different than South Bering (p = 0.00769)
## Cook Inlet is significantly different than South Bering (p = 0.0435)

# Add x- and y-coordinates for the comparison brackets
shan_dunn_plot_locale <- shan_dunn_results_locale %>%
  add_xy_position(x = "Locale")

# Inspect the table used for plotting
shan_dunn_plot_locale

# Box-and-whisker plot, with individual samples overlaid
shan_plot_locale <- ggplot(
  metadata,
  aes(x = Locale, y = Shannon, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_pvalue_manual(
    shan_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Locale",
    y = "Shannon Index per sample",
    title = "Shannon Index by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

shan_plot_locale

# gets both global and by-host summary data
summary(metadata$Shannon)

metadata %>%
  group_by(Locale) %>%
  summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_shannon))

### biological interpretation and conclusion
# Shannon entropy is significantly different in the Arctic and Cook Inlet compared to the South Bering. 
## This trend may be mostly driven by differences in Shannon entropy in beluga whales. 
## by host should be examined. 

# BY AGE---------------------------------------------------------------------------------
shan_test_age <- kruskal.test(Shannon ~ Age_Group, data = metadata_age)
shan_test_age #  p-value = 0.9588 not significant

# Dunn test to determine which groups are significantly different
shan_dunn_results_age<- metadata %>%
  dunn_test(
    Shannon ~ Age_Group,
    p.adjust.method = "hochberg"
  )
shan_dunn_results_age

### SIGNIFICANTLY DIFFERENT SHANNON ENTROPY
## there is no significant difference in shannon entropy between pups and non-pips
## in ringed and bearded seals

# Add x- and y-coordinates for the comparison brackets
shan_dunn_plot_age <- shan_dunn_results_age %>%
  add_xy_position(x = "Age_Group")

# Inspect the table used for plotting
shan_dunn_plot_age

# Box-and-whisker plot, with individual samples overlaid
shan_plot_age<- ggplot(
  metadata_age,
  aes(x = Age_Group, y = Shannon, fill = Predator)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  stat_pvalue_manual(
    shan_dunn_plot_age,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Age_Group",
    y = "Shannon Index per sample",
    title = "Shannon Index by Age_Group"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

shan_plot_age

# gets both global and by-host summary data
summary(metadata_age$Shannon)

metadata %>%
  group_by(Age_Group) %>%
  summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_shannon))

### biological interpretation and conclusion
# There are no significant differences in Shannon entropy 
## though i think there could be if Host is separated out, especially amongst ringed seals

# SHANNON PLOTS --------------------------------------------------------------

shan_plots <- (shan_plot_host + shan_plot_age + shan_plot_locale) / shan_plot_season / shan_plot_location  &
  theme(
    text = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

ggsave("Deliverables/ubiome/alphadiv/SHAN_ubiome-majorclass.png", plot = shan_plots, width = 40, height = 40, units = "in", dpi = 300)


# interpretation of all metrics collectively
### NEED TO CONTROL FOR HOST FOR ALL COVARIATES, OR RUN NEW TESTS. 

# HOST
host_plots <- obs_plot_host | bp_plot_host / shan_plot_host
host_plots

# SEASON
season_plots <- obs_plot_season | bp_plot_season / shan_plot_season
season_plots

# LOCATION
location_plots <- obs_plot_location | bp_plot_location / shan_plot_location
location_plots

# LOCALE
locale_plots <- obs_plot_locale | bp_plot_locale / shan_plot_locale
locale_plots

# AGE
age_plots <- obs_plot_age | bp_plot_age / shan_plot_age
age_plots


