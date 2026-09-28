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
load("./Scripts/12s/rdata/ps.12s.dietdat.Rdata")
load("./Scripts/ubiome/rdata/ubiome.ps.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.norep), "matrix")

# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------

# formatting for metadata
metadata <- data.frame(sample_data(ps.norep))

class(metadata)

# optionally removes fin whale, since there is only one sample (two sub-samples of the same sample)
metadata <- metadata %>%
  dplyr::filter(Predator != "fin whale")

asv_mat <- as.data.frame(asv_mat) %>% 
  rownames_to_column(var = "row_id") %>% 
  dplyr::filter(!row_id %in% c("WADE-003-150-A", "WADE-003-150-B")) %>% 
  column_to_rownames(var = "row_id")


# optionally removes spring and winter, since there are few samples
metadata <- metadata %>%
  dplyr::filter(season != "Winter")%>%
  dplyr::filter(season != "Spring")
asv_mat <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% metadata$LabID) %>%
  column_to_rownames(var = "LabID")

# optionally removes south bering from ringed seals (only 2 samples)
metadata <- metadata %>%
  dplyr::filter(
    !(Predator == "ringed seal" & Locale == "South Bering")
  )
asv_mat <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% metadata$LabID) %>%
  column_to_rownames(var = "LabID")


# formatting metadata for age class (only bearded and ringed seals)
metadata_age <- metadata %>%
  dplyr::filter(Predator %in% c("ringed seal", "bearded seal")) %>%
  dplyr::mutate(
    Age_Group = dplyr::case_when(
      Age_Group %in% c("", "Pending") ~ "Unknown",
      TRUE ~ Age_Group
    )
  )
# Optionally removed Unknown category: 
metadata_age <- metadata_age %>%
  dplyr::filter(!Age_Group == "Unknown")

# formatting asv_mat for age class (only bearded and ringed seals)
asv_mat_age <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

# sets host colors
host_colors <- c(
  "bearded seal" = "#F8766D",  # ggplot2 default red
  "beluga whale" = "#00BA38",  # ggplot2 default green
  "ringed seal"  = "#619CFF"   # ggplot2 default blue
)

# ------------------------------------------------------------------------------
# Gets richness, dominance, and information metrics for alpha diversity
# ------------------------------------------------------------------------------

alpha_div <- data.frame(
  LabID = rownames(asv_mat),
  
  # Observed ASV richness
  Observed = rowSums(asv_mat > 0),
  
  # Berger-Parker dominance: dominant ASV reads / total sample reads
  Berger_Parker = apply(asv_mat, 1, max) / rowSums(asv_mat),
  
  # Shannon diversity
  Shannon = vegan::diversity(asv_mat, index = "shannon"),
  
  # Useful QC covariate
  n_reads = rowSums(asv_mat),
  
  row.names = NULL
)

# adds alpha_div metrics to the metadata tables
metadata <- metadata %>%
  dplyr::left_join(alpha_div, by = "LabID")

metadata_age <- metadata_age %>%
  dplyr::left_join(alpha_div, by = "LabID")

# leftjoin diet alphadiv info to metadata

metadata_12s <- metadata_12s%>%
  dplyr::rename(diet_observed = Observed)%>%
  dplyr::rename(diet_bp = Berger_Parker)%>%
  dplyr::rename(diet_shan = Shannon)


metadata <- metadata %>%
  dplyr::left_join(
    metadata_12s %>%
      dplyr::select(LabID, diet_observed),
    by = "LabID"
  )

# ------------------------------------------------------------------------------
# Richness
# ------------------------------------------------------------------------------

# Robbins
## the ratio of singletons to total taxa


# Observed 
## the number of ASVs present

hist(metadata$Observed)

## wondering if I should run a Kruskal-Wallis test to test for actual differences here

# test for normalcy
norm <- shapiro.test(metadata$Observed)
norm # p-value = 1.755e-06 (significant - not normally distributed)

# test homogeneity of variance FOR HOST
#?leveneTest()
var <- leveneTest(Observed ~ Predator, metadata)
var # p-value = 4.613e-08 (significant - non-homonegeous variance)
# Does not meet assumptions for parametric tests

# tests for correlations between library depth and observed richness
cor.test(
  metadata$Observed,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
)

ggplot(metadata, aes(x = n_reads, y = Observed, color = Predator)) +
  geom_point(alpha = 0.75) +
  scale_x_log10() +
  theme_classic() +
  labs(
    x = "Retained reads per sample (log10 scale)",
    y = "Observed ASVs"
  )

# p-value = 0.5969; though there is a slightly positive trend, there is no evidence
# that the read depth is highly correlated with the observed richness
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)



# BY PREDATOR----------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_host <- kruskal.test(Observed ~ Predator, data = metadata)
obs_test_host #  p-value = 1.472e-08 (significant)

# Dunn test to determine which groups are significantly different
#?dunn_test()
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
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
## bearded seal and beluga whale (adjusted p-value = 0.0000000901)
## ringed seal and beluga whale (adjusted p-value = 0.0000141)

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

range(metadata$Year)

### biological interpretation and conclusion
# Bearded-seal and ringed-seal samples had higher observed ASV richness 
# than beluga-whale samples. There was no evidence of a difference in observed richness 
# between the two seal species

# templates for better wording of results!!!!
# Bacterial-community diversity and composition differed primarily among host 
# species and secondarily by season, with [host] exhibiting the highest alpha 
# diversity and [season] associated with [direction or community shift].

## for when I have R^2s for each to compare!
# Host species explained more variation in community composition than season 
# (r^2 vs r^2 vs r^2), while paired diet data identified [prey group or diet metric] as a 
# significant predictor of microbiome variation.


# ------------------------------------------------------------------------------
# Separates into Individual Host 
## because host is significant
# ------------------------------------------------------------------------------

beluga_metadata <- metadata %>%
  dplyr::filter(Predator == "beluga whale")
beluga_asv_mat <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% beluga_metadata$LabID) %>%
  column_to_rownames(var = "LabID")


ringed_metadata <- metadata %>%
  dplyr::filter(Predator == "ringed seal")
ringed_asv_mat <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% ringed_metadata$LabID) %>%
  column_to_rownames(var = "LabID")

bearded_metadata <- metadata %>%
  dplyr::filter(Predator == "bearded seal")
bearded_asv_mat <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% bearded_metadata$LabID) %>%
  column_to_rownames(var = "LabID")

# For the age group
ringed_metadata_age <- metadata_age %>%
  dplyr::filter(Predator == "ringed seal")
ringed_asv_mat_age <- asv_mat_age %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% ringed_metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

bearded_metadata_age <- metadata_age %>%
  dplyr::filter(Predator == "bearded seal")
bearded_asv_mat_age <- asv_mat_age %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% bearded_metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

# BY 12S DIET--------------------------------------------------------------------------------

# filters metadata to exclude NAs
beluga_metadata_diet <- metadata %>%
  dplyr::filter(!is.na(diet_observed))
beluga_asv_mat_diet <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% beluga_metadata_diet$LabID) %>%
  column_to_rownames(var = "LabID")

bearded_metadata_diet <- metadata %>%
  dplyr::filter(!is.na(diet_observed))
bearded_asv_mat_diet <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% bearded_metadata_diet$LabID) %>%
  column_to_rownames(var = "LabID")

ringed_metadata_diet <- metadata %>%
  dplyr::filter(!is.na(diet_observed))

ringed_asv_mat_diet <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% ringed_metadata_diet$LabID) %>%
  column_to_rownames(var = "LabID")


# Kruskall-Wallis test

# BELUGA
beluga_obs_test_diet <- kruskal.test(Observed ~ diet_observed, data = beluga_metadata_diet)
beluga_obs_test_diet # p-value = p-value = 0.05997 not significant 

# Dunn test to determine which groups are significantly different
beluga_obs_dunn_results_diet <- beluga_metadata_diet %>%
  dunn_test(
    Observed ~ diet_observed,
    p.adjust.method = "hochberg"
  )
beluga_obs_dunn_results_diet

# Add x- and y-coordinates for the comparison brackets
beluga_obs_dunn_plot_diet <- beluga_obs_dunn_results_diet %>%
  add_xy_position(x = "diet")

# Inspect the table used for plotting
beluga_obs_dunn_plot_diet


# filters metadata to exclude NAs
beluga_metadata_diet <- metadata %>%
  dplyr::filter(!is.na(diet_observed))

beluga_asv_mat_diet <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% beluga_metadata_diet$LabID) %>%
  column_to_rownames(var = "LabID")
# Kruskall-Wallis test

# bearded seal
bearded_obs_test_diet <- kruskal.test(Observed ~ diet_observed, data = bearded_metadata_diet)
bearded_obs_test_diet # p-value = p-value = 0.05997 not significant 

# Dunn test to determine which groups are significantly different
bearded_obs_dunn_results_diet <- bearded_metadata_diet %>%
  dunn_test(
    Observed ~ diet_observed,
    p.adjust.method = "hochberg"
  )
bearded_obs_dunn_results_diet

# Add x- and y-coordinates for the comparison brackets
bearded_obs_dunn_plot_diet <- bearded_obs_dunn_results_diet %>%
  add_xy_position(x = "diet")

# Inspect the table used for plotting
bearded_obs_dunn_plot_diet

# ringed seal
ringed_obs_test_diet <- kruskal.test(Observed ~ diet_observed, data = ringed_metadata_diet)
ringed_obs_test_diet # p-value = p-value = 0.05997 not significant 

# Dunn test to determine which groups are significantly different
ringed_obs_dunn_results_diet <- ringed_metadata_diet %>%
  dunn_test(
    Observed ~ diet_observed,
    p.adjust.method = "hochberg"
  )
ringed_obs_dunn_results_diet

# Add x- and y-coordinates for the comparison brackets
ringed_obs_dunn_plot_diet <- ringed_obs_dunn_results_diet %>%
  add_xy_position(x = "diet")

# Inspect the table used for plotting
ringed_obs_dunn_plot_diet


# BY SEASON---------------------------------------------------------------------------------
# Kruskall-Wallis test

# BELUGA
beluga_obs_test_season <- kruskal.test(Observed ~ season, data = beluga_metadata)
beluga_obs_test_season # p-value = p-value = 0.1524 not significant 

# Dunn test to determine which groups are significantly different
beluga_obs_dunn_results_season <- beluga_metadata %>%
  dunn_test(
    Observed ~ season,
    p.adjust.method = "hochberg"
  )
beluga_obs_dunn_results_season

# Add x- and y-coordinates for the comparison brackets
beluga_obs_dunn_plot_season <- beluga_obs_dunn_results_season %>%
  add_xy_position(x = "season")

# Inspect the table used for plotting
beluga_obs_dunn_plot_season


# RINGED
ringed_obs_test_season <- kruskal.test(Observed ~ season, data = ringed_metadata)
ringed_obs_test_season # p-value = p-value = 0.3935 not significant 

# Dunn test to determine which groups are significantly different
ringed_obs_dunn_results_season <- ringed_metadata %>%
  dunn_test(
    Observed ~ season,
    p.adjust.method = "hochberg"
  )
ringed_obs_dunn_results_season

# Add x- and y-coordinates for the comparison brackets
ringed_obs_dunn_plot_season <- ringed_obs_dunn_results_season %>%
  add_xy_position(x = "season")

# Inspect the table used for plotting
ringed_obs_dunn_plot_season

# BEARDED
bearded_obs_test_season <- kruskal.test(Observed ~ season, data = bearded_metadata)
bearded_obs_test_season # p-value = p-value = 0.7952 not significant 

# Dunn test to determine which groups are significantly different
bearded_obs_dunn_results_season <- bearded_metadata %>%
  dunn_test(
    Observed ~ season,
    p.adjust.method = "hochberg"
  )
bearded_obs_dunn_results_season

# Add x- and y-coordinates for the comparison brackets
bearded_obs_dunn_plot_season <- bearded_obs_dunn_results_season %>%
  add_xy_position(x = "season")

# Inspect the table used for plotting
bearded_obs_dunn_plot_season

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY SEASON
## no significantly differences in observed richness by season

# Box-and-whisker plot, with individual samples overlaid
# BELUGA BY SEASON
beluga_obs_plot_season <- ggplot(
  beluga_metadata,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    beluga_obs_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Number of ASVs per beluga sample",
    title = "Beluga observed ASV richness by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "right")

beluga_obs_plot_season

# RINGED BY SEASON
ringed_obs_plot_season <- ggplot(
  ringed_metadata,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    ringed_obs_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Number of ASVs per ringed seal sample",
    title = "Ringed seal observed ASV richness by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "right")

ringed_obs_plot_season

# BEARDED BY SEASON
bearded_obs_plot_season <- ggplot(
  bearded_metadata,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    bearded_obs_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Number of ASVs per bearded seal sample",
    title = "Bearded seal observed ASV richness by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "right")

bearded_obs_plot_season

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


# BY LOCALE-----------------------------------------------------------------------
# Kruskal-Wallis test

# BELUGA
beluga_obs_test_locale <- kruskal.test(Observed ~ Locale, data = beluga_metadata)
beluga_obs_test_locale # p-value = 0.7504 not significant

# Dunn test to determine which groups are significantly different
beluga_obs_dunn_results_locale <- beluga_metadata %>%
  dunn_test(
    Observed ~ Locale,
    p.adjust.method = "hochberg"
  )
beluga_obs_dunn_results_locale

beluga_obs_dunn_plot_locale <- beluga_obs_dunn_results_locale %>%
  add_xy_position(x = "Locale")

# Inspect the table used for plotting
beluga_obs_dunn_plot_locale


# Kruskal-Wallis test
# RINGED
n_ringed_locales <- dplyr::n_distinct(ringed_metadata_locale$Locale)

if (n_ringed_locales >= 2) {
  
  ringed_obs_test_locale <- kruskal.test(
    Observed ~ Locale,
    data = ringed_metadata_locale
  )
  
  ringed_obs_test_locale
  
  ringed_obs_dunn_results_locale <- ringed_metadata_locale %>%
    rstatix::dunn_test(
      Observed ~ Locale,
      p.adjust.method = "hochberg"
    )
  
  ringed_obs_dunn_results_locale
  
  ringed_obs_dunn_plot_locale <- ringed_obs_dunn_results_locale %>%
    rstatix::add_xy_position(x = "Locale")
  
  ringed_obs_plot_locale <- ggplot(
    ringed_metadata_locale,
    aes(x = Locale, y = Observed, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    ggpubr::stat_pvalue_manual(
      ringed_obs_dunn_plot_locale,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Locale",
      y = "Observed richness per ringed seal sample",
      title = "Ringed seal observed richness by locale"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  bearded_obs_plot_locale
  
} else {
  
  message(
    "Ringed seal locale test skipped: only one Locale is represented ",
    "after dropping missing values."
  )
  
  ringed_obs_plot_locale <- ggplot(
    ringed_metadata_locale,
    aes(x = Locale, y = Observed, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Locale",
      y = "Obersved richness per ringed seal sample",
      title = "Ringed seal observed richness by locale"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

ringed_obs_plot_locale

# BEARDED
#bearded_obs_test_locale <- kruskal.test(Observed ~ Locale, data = bearded_metadata)
## CAN'T RUN: ALL BEARDED SEAL SAMPLES ARE FROM THE ARCTIC DB
n_bearded_locales <- dplyr::n_distinct(bearded_metadata_locale$Locale)

if (n_bearded_locales >= 2) {
  
  bearded_obs_test_locale <- kruskal.test(
    Observed ~ Locale,
    data = bearded_metadata_locale
  )
  
  bearded_obs_test_locale
  
  bearded_obs_dunn_results_locale <- bearded_metadata_locale %>%
    rstatix::dunn_test(
      Observed ~ Locale,
      p.adjust.method = "hochberg"
    )
  
  bearded_obs_dunn_results_locale
  
  bearded_obs_dunn_plot_locale <- bearded_obs_dunn_results_locale %>%
    rstatix::add_xy_position(x = "Locale")
  
  bearded_obs_plot_locale <- ggplot(
    bearded_metadata_locale,
    aes(x = Locale, y = Observed, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    ggpubr::stat_pvalue_manual(
      bearded_obs_dunn_plot_locale,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Locale",
      y = "Observed richness per bearded seal sample",
      title = "Bearded seal observed richness by locale"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  bearded_obs_plot_locale
  
} else {
  
  message(
    "Bearded seal locale test skipped: only one Locale is represented ",
    "after dropping missing values."
  )
  
  bearded_obs_plot_locale <- ggplot(
    bearded_metadata_locale,
    aes(x = Locale, y = Observed, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Locale",
      y = "Obersved richness per bearded seal sample",
      title = "Bearded seal observed richness by locale"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

bearded_obs_plot_locale

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY LOCALE
# NONE

# Box-and-whisker plot, with individual samples overlaid
# BY LOCALE
beluga_obs_plot_locale <- ggplot(
  beluga_metadata,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    beluga_obs_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Locale",
    y = "Number of ASVs per beluga sample",
    title = "Beluga observed ASV richness by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

beluga_obs_plot_locale

ringed_obs_plot_locale <- ggplot(
  ringed_metadata,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    ringed_obs_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Locale",
    y = "Number of ASVs per ringed sample",
    title = "ringed seal observed ASV richness by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ringed_obs_plot_locale

bearded_obs_plot_locale <- ggplot(
  bearded_metadata,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  # stat_pvalue_manual(
  #   bearded_obs_dunn_plot_locale,
  #   label = "p.adj.signif",
  #   hide.ns = TRUE,
  #   tip.length = 0.01
  # ) +
  labs(
    x = "Locale",
    y = "Number of ASVs per ringed sample",
    title = "bearded seal observed ASV richness by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

bearded_obs_plot_locale

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
# NONE

# BY AGE--------------------------------------------------------------------------
# Kruskal-Wallis test
# RINGED
ringed_obs_test_age <- kruskal.test(Observed ~ Age_Group, data = ringed_metadata_age)
ringed_obs_test_age # p-value = 0.4567 not significant

# Dunn test to determine which groups are significantly different
ringed_obs_dunn_results_age <- ringed_metadata_age %>%
  dunn_test(
    Observed ~ Age_Group,
    p.adjust.method = "hochberg"
  )
ringed_obs_dunn_results_age

ringed_obs_dunn_plot_age <- ringed_obs_dunn_results_age %>%
  add_xy_position(x = "Age_Group")

# Inspect the table used for plotting
ringed_obs_dunn_plot_age

# BEARDED
bearded_obs_test_age <- kruskal.test(Observed ~ Age_Group, data = bearded_metadata_age)
bearded_obs_test_age # p-value = 0.4249 not significant

# Dunn test to determine which groups are significantly different
bearded_obs_dunn_results_age <- bearded_metadata_age %>%
  dunn_test(
    Observed ~ Age_Group,
    p.adjust.method = "hochberg"
  )
bearded_obs_dunn_results_age

bearded_obs_dunn_plot_age <- bearded_obs_dunn_results_age %>%
  add_xy_position(x = "Age_Group")

# Inspect the table used for plotting
bearded_obs_dunn_plot_age

### SIGNIFICANTLY DIFFERENT OBSERVED RICHNESS BY AGE
# none

# Box-and-whisker plot, with individual samples overlaid
# BY AGE GROUP
ringed_obs_plot_age <- ggplot(
  ringed_metadata_age,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    ringed_obs_dunn_plot_age,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Age_Group",
    y = "Number of ASVs per ringed seal sample",
    title = "Ringed seal observed ASV richness by Age Group"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ringed_obs_plot_age

bearded_obs_plot_age <- ggplot(
  bearded_metadata_age,
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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  stat_pvalue_manual(
    bearded_obs_dunn_plot_age,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Age_Group",
    y = "Number of ASVs per ringed seal sample",
    title = "Beardeds seal observed ASV richness by Age Group"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

bearded_obs_plot_age

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

# Host-level dominance plot
obs_plot_host

# Seasonal dominance plots, one per host
obs_plots_season <- (
  bearded_obs_plot_season+
  beluga_obs_plot_season +
  ringed_obs_plot_season 
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

obs_plots_season

# Locale obs plots, one per host
obs_plots_locale <- (
  bearded_obs_plot_locale+
  beluga_obs_plot_locale +
    ringed_obs_plot_locale
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

obs_plots_locale


# Known-age obs plots for ringed and bearded seals
obs_plots_age <- (
  bearded_obs_plot_age+
  ringed_obs_plot_age
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

obs_plots_age


# Combined dominance figure
obs_plots <- (
  obs_plot_host /
    obs_plots_season /
    # obs_plots_locale /
    obs_plots_age
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

obs_plots

ggsave(
  "Deliverables/ubiome/alphadiv/OBS_ubiome.png",
  plot = obs_plots,
  width = 20,
  height = 26,
  units = "in",
  dpi = 300
)

ggsave(
  "Deliverables/ubiome/alphadiv/OBS-HOST_ubiome.png",
  plot = obs_plot_host,
  width = 20,
  height = 26,
  units = "in",
  dpi = 300
)

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
bp.norm # p-value = 1.369e-07 (significant; nonnormal)

# test homogeneity of variance
#?leveneTest()
bp.var <- leveneTest(Berger_Parker ~ Predator, metadata)
bp.var # p-value = 0.8509 (non-signficant - homogeneous variance)

# Does meet assumptions for parametric tests

# tests if read depth is correlated with Berger-Parker ratios
cor.test(
  metadata$Berger_Parker,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
) 

# p-value = 0.8443; though there is a slight negative trend, there is no evidence
# that the read depth is highly correlated with the Berger-Parker ratio
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)

# BY HOST---------------------------------------------------------------------------------
# Kruskall-Wallis test
bp_test_host <- kruskal.test(Berger_Parker ~ Predator, data = metadata)
bp_test_host #  p-value = 0.005746 significant 

# Dunn test to determine which groups are significantly different
bp_dunn_results_host <- metadata %>%
  dunn_test(
    Berger_Parker ~ Predator,
    p.adjust.method = "hochberg"
  )
bp_dunn_results_host

### SIGNIFICANTLY DIFFERENT DOMINANCE
## bearded seal and beluga whale (adjusted p-value = 0.00495)

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
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
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

# SEASON-----------------------------------------------------------------------------------------
# Remove samples with missing, blank, or whitespace-only season data separately
# within each host.

beluga_metadata_season <- beluga_metadata %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    season = droplevels(factor(season))
  )

ringed_metadata_season <- ringed_metadata %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    season = droplevels(factor(season))
  )

bearded_metadata_season <- bearded_metadata %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    season = droplevels(factor(season))
  )

# Check seasonal sample sizes before interpreting tests
beluga_metadata_season %>%
  dplyr::count(season, name = "n")

ringed_metadata_season %>%
  dplyr::count(season, name = "n")

bearded_metadata_season %>%
  dplyr::count(season, name = "n")

# BELUGA: Berger-Parker dominance by season
beluga_bp_test_season <- kruskal.test(
  Berger_Parker ~ season,
  data = beluga_metadata_season
)

beluga_bp_test_season # p-val ue = 0.01515 not signficant

beluga_bp_dunn_results_season <- beluga_metadata_season %>%
  rstatix::dunn_test(
    Berger_Parker ~ season,
    p.adjust.method = "hochberg"
  )

beluga_bp_dunn_results_season

beluga_bp_dunn_plot_season <- beluga_bp_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

beluga_bp_dunn_plot_season

beluga_bp_plot_season <- ggplot(
  beluga_metadata_season,
  aes(x = season, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  ggpubr::stat_pvalue_manual(
    beluga_bp_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Dominance of ASVs per beluga sample",
    title = "Beluga ASV dominance by season"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

beluga_bp_plot_season

# RINGED SEAL: Berger-Parker dominance by season
ringed_bp_test_season <- kruskal.test(
  Berger_Parker ~ season,
  data = ringed_metadata_season
)

ringed_bp_test_season # p-value = 0.349 not significant

ringed_bp_dunn_results_season <- ringed_metadata_season %>%
  rstatix::dunn_test(
    Berger_Parker ~ season,
    p.adjust.method = "hochberg"
  )

ringed_bp_dunn_results_season

ringed_bp_dunn_plot_season <- ringed_bp_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

ringed_bp_dunn_plot_season

ringed_bp_plot_season <- ggplot(
  ringed_metadata_season,
  aes(x = season, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  ggpubr::stat_pvalue_manual(
    ringed_bp_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Dominance of ASVs per ringed seal sample",
    title = "Ringed seal ASV dominance by season"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ringed_bp_plot_season

# BEARDED SEAL: Berger-Parker dominance by season
bearded_bp_test_season <- kruskal.test(
  Berger_Parker ~ season,
  data = bearded_metadata_season
)

bearded_bp_test_season # p-value = 0.7066 not significant

bearded_bp_dunn_results_season <- bearded_metadata_season %>%
  rstatix::dunn_test(
    Berger_Parker ~ season,
    p.adjust.method = "hochberg"
  )

bearded_bp_dunn_results_season

bearded_bp_dunn_plot_season <- bearded_bp_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

bearded_bp_dunn_plot_season

bearded_bp_plot_season <- ggplot(
  bearded_metadata_season,
  aes(x = season, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  ggpubr::stat_pvalue_manual(
    bearded_bp_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Season",
    y = "Dominance of ASVs per bearded seal sample",
    title = "Bearded seal ASV dominance by season"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

bearded_bp_plot_season

# Summaries by host and season
bp_summary_season <- dplyr::bind_rows(
  beluga_metadata_season,
  ringed_metadata_season,
  bearded_metadata_season
) %>%
  dplyr::group_by(Predator, season) %>%
  dplyr::summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Predator, dplyr::desc(median_Berger_Parker))

bp_summary_season

# ------------------------------------------------------------------------------
# BY LOCALE
# ------------------------------------------------------------------------------

# Remove samples with missing, blank, or whitespace-only Locale data separately
# within each host.

beluga_metadata_locale <- beluga_metadata %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    Locale = droplevels(factor(Locale))
  )

ringed_metadata_locale <- ringed_metadata %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    Locale = droplevels(factor(Locale))
  )

bearded_metadata_locale <- bearded_metadata %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    Locale = droplevels(factor(Locale))
  )

# Check locale sample sizes before testing
beluga_metadata_locale %>%
  dplyr::count(Locale, name = "n")

ringed_metadata_locale %>%
  dplyr::count(Locale, name = "n")

bearded_metadata_locale %>%
  dplyr::count(Locale, name = "n")

# BELUGA: Berger-Parker dominance by locale
beluga_bp_test_locale <- kruskal.test(
  Berger_Parker ~ Locale,
  data = beluga_metadata_locale
)

beluga_bp_test_locale # p-value = 0.01614 significant

beluga_bp_dunn_results_locale <- beluga_metadata_locale %>%
  rstatix::dunn_test(
    Berger_Parker ~ Locale,
    p.adjust.method = "hochberg"
  )

beluga_bp_dunn_results_locale 
## Artic and South Bering (0.0356)
## Cook Inlet and South Bering (0.0356)

beluga_bp_dunn_plot_locale <- beluga_bp_dunn_results_locale %>%
  rstatix::add_xy_position(x = "Locale")

beluga_bp_dunn_plot_locale

beluga_bp_plot_locale <- ggplot(
  beluga_metadata_locale,
  aes(x = Locale, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  ggpubr::stat_pvalue_manual(
    beluga_bp_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  labs(
    x = "Locale",
    y = "Dominance of ASVs per beluga sample",
    title = "Beluga ASV dominance by locale"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

beluga_bp_plot_locale

## biological interpretation
## beluga whales from the Arctic are much more dominated by one class (?) 
## of microbes than the South Bering is. The South Bering sea belugas have much 
## lower dominance/higher community evenness than the Arctic

# RINGED SEAL: Berger-Parker dominance by locale
# ringed_bp_test_locale <- kruskal.test(
#   Berger_Parker ~ Locale,
#   data = ringed_metadata_locale
# )
# 
# ringed_bp_test_locale # p-value = 0.6718 not significant
# 
# ringed_bp_dunn_results_locale <- ringed_metadata_locale %>%
#   rstatix::dunn_test(
#     Berger_Parker ~ Locale,
#     p.adjust.method = "hochberg"
#   )
# 
# ringed_bp_dunn_results_locale
# 
# ringed_bp_dunn_plot_locale <- ringed_bp_dunn_results_locale %>%
#   rstatix::add_xy_position(x = "Locale")
# 
# ringed_bp_dunn_plot_locale

ringed_bp_plot_locale <- ggplot(
  ringed_metadata_locale,
  aes(x = Locale, y = Berger_Parker, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  )+
  # ggpubr::stat_pvalue_manual(
  #   ringed_bp_dunn_plot_locale,
  #   label = "p.adj.signif",
  #   hide.ns = TRUE,
  #   tip.length = 0.01
  # ) +
  labs(
    x = "Locale",
    y = "Dominance of ASVs per ringed seal sample",
    title = "Ringed seal ASV dominance by locale"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ringed_bp_plot_locale

# BEARDED SEAL: Berger-Parker dominance by locale
#
# This host may have only one locale. Kruskal-Wallis and Dunn tests require at
# least two represented locale groups, so first check the number of levels.

n_bearded_locales <- dplyr::n_distinct(bearded_metadata_locale$Locale)

if (n_bearded_locales >= 2) {
  
  bearded_bp_test_locale <- kruskal.test(
    Berger_Parker ~ Locale,
    data = bearded_metadata_locale
  )
  
  bearded_bp_test_locale
  
  bearded_bp_dunn_results_locale <- bearded_metadata_locale %>%
    rstatix::dunn_test(
      Berger_Parker ~ Locale,
      p.adjust.method = "hochberg"
    )
  
  bearded_bp_dunn_results_locale
  
  bearded_bp_dunn_plot_locale <- bearded_bp_dunn_results_locale %>%
    rstatix::add_xy_position(x = "Locale")
  
  bearded_bp_plot_locale <- ggplot(
    bearded_metadata_locale,
    aes(x = Locale, y = Berger_Parker, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    )+
    ggpubr::stat_pvalue_manual(
      bearded_bp_dunn_plot_locale,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    labs(
      x = "Locale",
      y = "Dominance of ASVs per bearded seal sample",
      title = "Bearded seal ASV dominance by locale"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  bearded_bp_plot_locale
  
} else {
  
  message(
    "Bearded seal locale test skipped: only one Locale is represented ",
    "after dropping missing values."
  )
  
  bearded_bp_plot_locale <- ggplot(
    bearded_metadata_locale,
    aes(x = Locale, y = Berger_Parker, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    )+
    labs(
      x = "Locale",
      y = "Dominance of ASVs per bearded seal sample",
      title = "Bearded seal ASV dominance by locale"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

bearded_bp_plot_locale

# Summaries by host and locale
bp_summary_locale <- dplyr::bind_rows(
  beluga_metadata_locale,
  ringed_metadata_locale,
  bearded_metadata_locale
) %>%
  dplyr::group_by(Predator, Locale) %>%
  dplyr::summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Predator, dplyr::desc(median_Berger_Parker))

bp_summary_locale

# ------------------------------------------------------------------------------
# BY AGE GROUP
# ------------------------------------------------------------------------------

# Age analysis uses ringed and bearded seals only.
# Unknown age is excluded because it is missing information, not an age class.
# NA, blank, Pending, and Unknown values are not included in the primary test.

ringed_metadata_age_known <- ringed_metadata_age %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(Age_Group),
    trimws(as.character(Age_Group)) != "",
    Age_Group %in% c("Non-pup", "Pup")
  ) %>%
  dplyr::mutate(
    Age_Group = droplevels(
      factor(Age_Group, levels = c("Non-pup", "Pup"))
    )
  )

bearded_metadata_age_known <- bearded_metadata_age %>%
  dplyr::filter(
    !is.na(Berger_Parker),
    !is.na(Age_Group),
    trimws(as.character(Age_Group)) != "",
    Age_Group %in% c("Non-pup", "Pup")
  ) %>%
  dplyr::mutate(
    Age_Group = droplevels(
      factor(Age_Group, levels = c("Non-pup", "Pup"))
    )
  )

# Check known-age sample sizes
ringed_metadata_age_known %>%
  dplyr::count(Age_Group, name = "n")

bearded_metadata_age_known %>%
  dplyr::count(Age_Group, name = "n")

# RINGED SEAL: Berger-Parker dominance by age group
n_ringed_age_groups <- dplyr::n_distinct(ringed_metadata_age_known$Age_Group)

  ringed_bp_test_age <- kruskal.test(
    Berger_Parker ~ Age_Group,
    data = ringed_metadata_age_known
  )
  
  ringed_bp_test_age # p-value = 0.3673 not significant
  
  ringed_bp_dunn_results_age <- ringed_metadata_age_known %>%
    rstatix::dunn_test(
      Berger_Parker ~ Age_Group,
      p.adjust.method = "hochberg"
    )
  
  ringed_bp_dunn_results_age
  
  ringed_bp_dunn_plot_age <- ringed_bp_dunn_results_age %>%
    rstatix::add_xy_position(x = "Age_Group")
  
  ringed_bp_plot_age <- ggplot(
    ringed_metadata_age_known,
    aes(x = Age_Group, y = Berger_Parker, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    )+
    ggpubr::stat_pvalue_manual(
      ringed_bp_dunn_plot_age,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    labs(
      x = "Age group",
      y = "Dominance of ASVs per ringed seal sample",
      title = "Ringed seal ASV dominance by age group"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  ringed_bp_plot_age

  
  ringed_bp_plot_age <- ggplot(
    ringed_metadata_age_known,
    aes(x = Age_Group, y = Berger_Parker, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    )+
    labs(
      x = "Age group",
      y = "Dominance of ASVs per ringed seal sample",
      title = "Ringed seal ASV dominance by age group"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")

ringed_bp_plot_age

# BEARDED SEAL: Berger-Parker dominance by age group
n_bearded_age_groups <- dplyr::n_distinct(bearded_metadata_age_known$Age_Group)

  bearded_bp_test_age <- kruskal.test(
    Berger_Parker ~ Age_Group,
    data = bearded_metadata_age_known
  )
  
  bearded_bp_test_age # p-value = 0.3673 not significant
  
  bearded_bp_dunn_results_age <- bearded_metadata_age_known %>%
    rstatix::dunn_test(
      Berger_Parker ~ Age_Group,
      p.adjust.method = "hochberg"
    )
  
  bearded_bp_dunn_results_age
  
  bearded_bp_dunn_plot_age <- bearded_bp_dunn_results_age %>%
    rstatix::add_xy_position(x = "Age_Group")
  
  bearded_bp_plot_age <- ggplot(
    bearded_metadata_age_known,
    aes(x = Age_Group, y = Berger_Parker, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    )+
    ggpubr::stat_pvalue_manual(
      bearded_bp_dunn_plot_age,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    labs(
      x = "Age group",
      y = "Dominance of ASVs per bearded seal sample",
      title = "Bearded seal ASV dominance by age group"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  bearded_bp_plot_age
  
  
  bearded_bp_plot_age <- ggplot(
    bearded_metadata_age_known,
    aes(x = Age_Group, y = Berger_Parker, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    )+
    labs(
      x = "Age group",
      y = "Dominance of ASVs per bearded seal sample",
      title = "Bearded seal ASV dominance by age group"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")

bearded_bp_plot_age

# Summaries by host and known age group
bp_summary_age <- dplyr::bind_rows(
  ringed_metadata_age_known,
  bearded_metadata_age_known
) %>%
  dplyr::group_by(Predator, Age_Group) %>%
  dplyr::summarise(
    n = sum(!is.na(Berger_Parker)),
    median_Berger_Parker = median(Berger_Parker, na.rm = TRUE),
    mean_Berger_Parker = mean(Berger_Parker, na.rm = TRUE),
    IQR_Berger_Parker = IQR(Berger_Parker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Predator, dplyr::desc(median_Berger_Parker))

bp_summary_age

# ------------------------------------------------------------------------------
# DOMINANCE PLOTS
# ------------------------------------------------------------------------------

# Host-level dominance plot
bp_plot_host

# Seasonal dominance plots, one per host
bp_plots_season <- (
  bearded_bp_plot_season+
  beluga_bp_plot_season +
    ringed_bp_plot_season
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

bp_plots_season

# Locale dominance plots, one per host
bp_plots_locale <- (
  bearded_bp_plot_locale+
  beluga_bp_plot_locale +
    ringed_bp_plot_locale
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

bp_plots_locale


# Known-age dominance plots for ringed and bearded seals
bp_plots_age <- (
  bearded_bp_plot_age+
  ringed_bp_plot_age
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

bp_plots_age


# Combined dominance figure
bp_plots <- (
  bp_plot_host /
    bp_plots_season /
    bp_plots_locale /
    bp_plots_age
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

bp_plots

ggsave(
  "Deliverables/ubiome/alphadiv/BPDOM_ubiome-majorclass.png",
  plot = bp_plots,
  width = 20,
  height = 26,
  units = "in",
  dpi = 300
)

# DOMINANT ASV -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# GET TAXONOMY AND DOMINANT ASV FOR EACH SAMPLE
# ------------------------------------------------------------------------------

# Extract taxonomy table and retain ASV IDs as a column
tax_df <- as.data.frame(phyloseq::tax_table(ps.norep)) %>%
  tibble::rownames_to_column("ASV")

# Identify the column index of the most abundant ASV for each sample.
# ties.method = "first" consistently selects the first ASV if two ASVs tie.
dominant_idx <- max.col(asv_mat, ties.method = "first")

# Create a sample-level table containing each sample's dominant ASV.
dominant_asv_per_sample <- tibble::tibble(
  LabID = rownames(asv_mat),
  ASV = colnames(asv_mat)[dominant_idx],
  n_reads_asv = asv_mat[
    cbind(seq_len(nrow(asv_mat)), dominant_idx)
  ],
  sample_reads = rowSums(asv_mat)
) %>%
  dplyr::mutate(
    dominance = n_reads_asv / sample_reads
  ) %>%
  dplyr::left_join(
    tax_df,
    by = "ASV"
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(
        LabID,
        Predator,
        season,
        Locale,
        Age_Group
      ),
    by = "LabID"
  ) %>%
  dplyr::select(
    LabID,
    Predator,
    season,
    Locale,
    Age_Group,
    ASV,
    n_reads_asv,
    sample_reads,
    dominance,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus
  ) %>%
  dplyr::arrange(dplyr::desc(dominance))

dominant_asv_per_sample

# ------------------------------------------------------------------------------
# CREATE READABLE DOMINANT-ASV LABELS
# ------------------------------------------------------------------------------

# ASV sequences or ASV IDs are often difficult to read in a figure legend.
# This label uses genus where available, otherwise Family, Order, Class, or
# Phylum, and appends the ASV ID to distinguish multiple ASVs in the same taxon.

dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::mutate(
    ASV_short = substr(ASV, 1, 12),
    
    taxon_label = dplyr::case_when(
      !is.na(Genus) &
        Genus != "" &
        !tolower(Genus) %in% c("na", "unknown", "unclassified") ~
        paste0("g__ ", Genus, " | ", ASV_short),
      
      !is.na(Family) &
        Family != "" &
        !tolower(Family) %in% c("na", "unknown", "unclassified") ~
        paste0("f__ ", Family, " | ", ASV_short),
      
      !is.na(Order) &
        Order != "" &
        !tolower(Order) %in% c("na", "unknown", "unclassified") ~
        paste0("o__ ", Order, " | ", ASV_short),
      
      !is.na(Class) &
        Class != "" &
        !tolower(Class) %in% c("na", "unknown", "unclassified") ~
        paste0("c__ ", Class, " | ", ASV_short),
      
      !is.na(Phylum) &
        Phylum != "" &
        !tolower(Phylum) %in% c("na", "unknown", "unclassified") ~
        paste0("p__ ", Phylum, " | ", ASV_short),
      
      TRUE ~ paste0("Unclassified | ", ASV_short)
    )
  )

# ------------------------------------------------------------------------------
# OVERALL DOMINANT-ASV SUMMARY
# ------------------------------------------------------------------------------

dominant_asv_summary <- dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

dominant_asv_summary

write.csv(
  dominant_asv_summary,
  "Deliverables/ubiome/dominant_asv_summary_all_hosts.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# SPLIT DOMINANT-ASV DATA BY HOST
# ------------------------------------------------------------------------------

beluga_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(Predator == "beluga whale")

ringed_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(Predator == "ringed seal")

bearded_dominant_asv_per_sample <- dominant_asv_per_sample %>%
  dplyr::filter(Predator == "bearded seal")

# For age analysis, retain only seal samples and then create separate host data.
# The primary age analysis uses known Non-pup and Pup values only.
ringed_dominant_asv_per_sample_age <- dominant_asv_per_sample %>%
  dplyr::filter(
    Predator == "ringed seal",
    Age_Group %in% c("Non-pup", "Pup")
  )

bearded_dominant_asv_per_sample_age <- dominant_asv_per_sample %>%
  dplyr::filter(
    Predator == "bearded seal",
    Age_Group %in% c("Non-pup", "Pup")
  )

# Confirm that dominant-ASV sample tables match their host metadata objects.
setdiff(beluga_dominant_asv_per_sample$LabID, beluga_metadata$LabID)
setdiff(ringed_dominant_asv_per_sample$LabID, ringed_metadata$LabID)
setdiff(bearded_dominant_asv_per_sample$LabID, bearded_metadata$LabID)

# ------------------------------------------------------------------------------
# HOST-SPECIFIC DOMINANT-ASV SUMMARY TABLES
# ------------------------------------------------------------------------------

beluga_dominant_asv_summary <- beluga_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(beluga_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

ringed_dominant_asv_summary <- ringed_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(ringed_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

bearded_dominant_asv_summary <- bearded_dominant_asv_per_sample %>%
  dplyr::group_by(
    ASV,
    taxon_label,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus
  ) %>%
  dplyr::summarise(
    n_samples_dominant = dplyr::n(),
    proportion_samples_dominant = n_samples_dominant / nrow(bearded_dominant_asv_per_sample),
    mean_dominance = mean(dominance, na.rm = TRUE),
    median_dominance = median(dominance, na.rm = TRUE),
    max_dominance = max(dominance, na.rm = TRUE),
    total_reads_when_dominant = sum(n_reads_asv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_samples_dominant),
    dplyr::desc(mean_dominance)
  )

beluga_dominant_asv_summary
ringed_dominant_asv_summary
bearded_dominant_asv_summary

write.csv(
  beluga_dominant_asv_summary,
  "Deliverables/ubiome/beluga_dominant_asv_summary.csv",
  row.names = FALSE
)

write.csv(
  ringed_dominant_asv_summary,
  "Deliverables/ubiome/ringed_dominant_asv_summary.csv",
  row.names = FALSE
)

write.csv(
  bearded_dominant_asv_summary,
  "Deliverables/ubiome/bearded_dominant_asv_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# IDENTIFY TOP DOMINANT ASVs FOR EACH HOST
# ------------------------------------------------------------------------------

# To prevent unreadable stacked-bar legends, retain the top 8 dominant ASVs
# within each host. All remaining ASVs are grouped into "Other dominant ASVs".
# Change n = 8 if you want more or fewer individually displayed ASVs.

n_top_asvs <- 8

beluga_top_dominant_asvs <- beluga_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

ringed_top_dominant_asvs <- ringed_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

bearded_top_dominant_asvs <- bearded_dominant_asv_summary %>%
  dplyr::slice_head(n = n_top_asvs) %>%
  dplyr::pull(taxon_label)

# ------------------------------------------------------------------------------
# FUNCTION: PLOT DOMINANT-ASV IDENTITY WITHIN A HOST
# ------------------------------------------------------------------------------

# The plot displays:
#   x-axis: a grouping variable such as season, Locale, or Age_Group
#   y-axis: proportion of samples in that group
#   fill: ASV that had the highest within-sample read count
#
# position = "fill" makes every bar sum to 100%.

plot_dominant_asv <- function(
    data,
    grouping_var,
    host_name,
    top_asvs,
    x_label,
    title
) {
  
  grouping_var <- rlang::ensym(grouping_var)
  grouping_var_name <- rlang::as_string(grouping_var)
  
  plot_df <- data %>%
    dplyr::filter(
      !is.na(.data[[grouping_var_name]]),
      trimws(as.character(.data[[grouping_var_name]])) != ""
    ) %>%
    dplyr::mutate(
      dominant_asv_plot = dplyr::if_else(
        taxon_label %in% top_asvs,
        taxon_label,
        "Other dominant ASVs"
      )
    ) %>%
    dplyr::count(
      !!grouping_var,
      dominant_asv_plot,
      name = "n_samples"
    ) %>%
    dplyr::group_by(!!grouping_var) %>%
    dplyr::mutate(
      group_n = sum(n_samples),
      proportion_samples = n_samples / group_n
    ) %>%
    dplyr::ungroup()
  
  ggplot(
    plot_df,
    aes(
      x = !!grouping_var,
      y = proportion_samples,
      fill = dominant_asv_plot
    )
  ) +
    geom_col(
      position = "fill",
      width = 0.75,
      colour = "black",
      linewidth = 0.15
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_fill_manual(
      values = c(
        setNames(
          scales::hue_pal()(length(top_asvs)),
          top_asvs
        ),
        "Other dominant ASVs" = "grey70"
      ),
      drop = FALSE
    ) +
    labs(
      x = x_label,
      y = "Proportion of samples",
      fill = "Dominant ASV",
      title = title
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      legend.position = "right"
    )
}

# ------------------------------------------------------------------------------
# DOMINANT ASV BY SEASON, WITHIN EACH HOST
# ------------------------------------------------------------------------------

# BELUGA BY SEASON
beluga_dominant_asv_plot_season <- plot_dominant_asv(
  data = beluga_dominant_asv_per_sample,
  grouping_var = season,
  host_name = "beluga whale",
  top_asvs = beluga_top_dominant_asvs,
  x_label = "Season",
  title = "Dominant ASVs in beluga samples by season"
)

beluga_dominant_asv_plot_season

# RINGED SEAL BY SEASON
ringed_dominant_asv_plot_season <- plot_dominant_asv(
  data = ringed_dominant_asv_per_sample,
  grouping_var = season,
  host_name = "ringed seal",
  top_asvs = ringed_top_dominant_asvs,
  x_label = "Season",
  title = "Dominant ASVs in ringed seal samples by season"
)

ringed_dominant_asv_plot_season

# BEARDED SEAL BY SEASON
bearded_dominant_asv_plot_season <- plot_dominant_asv(
  data = bearded_dominant_asv_per_sample,
  grouping_var = season,
  host_name = "bearded seal",
  top_asvs = bearded_top_dominant_asvs,
  x_label = "Season",
  title = "Dominant ASVs in bearded seal samples by season"
)

bearded_dominant_asv_plot_season

# Combine seasonal dominant-ASV plots
dominant_asv_plots_season <- (
  bearded_dominant_asv_plot_season /
  beluga_dominant_asv_plot_season /
    ringed_dominant_asv_plot_season
    
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

dominant_asv_plots_season

ggsave(
  "Deliverables/ubiome/dominant_asv_by-season_each-host.png",
  plot = dominant_asv_plots_season,
  width = 18,
  height = 18,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# DOMINANT ASV BY LOCALE, WITHIN EACH HOST
# ------------------------------------------------------------------------------

# BELUGA BY LOCALE
beluga_dominant_asv_plot_locale <- plot_dominant_asv(
  data = beluga_dominant_asv_per_sample,
  grouping_var = Locale,
  host_name = "beluga whale",
  top_asvs = beluga_top_dominant_asvs,
  x_label = "Locale",
  title = "Dominant ASVs in beluga samples by locale"
)

beluga_dominant_asv_plot_locale

# RINGED SEAL BY LOCALE
ringed_dominant_asv_plot_locale <- plot_dominant_asv(
  data = ringed_dominant_asv_per_sample,
  grouping_var = Locale,
  host_name = "ringed seal",
  top_asvs = ringed_top_dominant_asvs,
  x_label = "Locale",
  title = "Dominant ASVs in ringed seal samples by locale"
)

ringed_dominant_asv_plot_locale

# BEARDED SEAL BY LOCALE
#
# All bearded seal samples may be Arctic. The plot remains descriptive, but a
# comparison among locales is not possible unless at least two Locale values exist.

bearded_dominant_asv_plot_locale <- plot_dominant_asv(
  data = bearded_dominant_asv_per_sample,
  grouping_var = Locale,
  host_name = "bearded seal",
  top_asvs = bearded_top_dominant_asvs,
  x_label = "Locale",
  title = "Dominant ASVs in bearded seal samples by locale"
)

bearded_dominant_asv_plot_locale

# Combine locale dominant-ASV plots
dominant_asv_plots_locale <- (
  bearded_dominant_asv_plot_locale / 
  beluga_dominant_asv_plot_locale /
    ringed_dominant_asv_plot_locale
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

dominant_asv_plots_locale

ggsave(
  "Deliverables/ubiome/dominant_asv_by-locale_each-host.png",
  plot = dominant_asv_plots_locale,
  width = 18,
  height = 18,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# DOMINANT ASV BY AGE GROUP, WITHIN EACH SEAL HOST
# ------------------------------------------------------------------------------

# RINGED SEAL BY AGE GROUP
ringed_dominant_asv_plot_age <- plot_dominant_asv(
  data = ringed_dominant_asv_per_sample_age,
  grouping_var = Age_Group,
  host_name = "ringed seal",
  top_asvs = ringed_top_dominant_asvs,
  x_label = "Age group",
  title = "Dominant ASVs in ringed seal samples by age group"
)

ringed_dominant_asv_plot_age

# BEARDED SEAL BY AGE GROUP
bearded_dominant_asv_plot_age <- plot_dominant_asv(
  data = bearded_dominant_asv_per_sample_age,
  grouping_var = Age_Group,
  host_name = "bearded seal",
  top_asvs = bearded_top_dominant_asvs,
  x_label = "Age group",
  title = "Dominant ASVs in bearded seal samples by age group"
)

bearded_dominant_asv_plot_age

# Combine age-group dominant-ASV plots
dominant_asv_plots_age <- (
  bearded_dominant_asv_plot_age / 
  ringed_dominant_asv_plot_age
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

dominant_asv_plots_age

ggsave(
  "Deliverables/ubiome/dominant_asv_by-age_each-seal-host.png",
  plot = dominant_asv_plots_age,
  width = 18,
  height = 9,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# HOST-SPECIFIC DOMINANT-ASV STRENGTH PLOTS
# ------------------------------------------------------------------------------

# These plots answer a different question than the stacked ASV identity plots:
# they show HOW strongly the dominant ASV dominates each sample.

beluga_dominance_plot_season <- ggplot(
  beluga_dominant_asv_per_sample %>%
    dplyr::filter(
      !is.na(season),
      trimws(as.character(season)) != ""
    ),
  aes(x = season, y = dominance, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    colour = "black"
  ) +
  geom_jitter(
    width = 0.15,
    shape = 21,
    colour = "black",
    size = 2,
    alpha = 0.8,
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Season",
    y = "Reads represented by dominant ASV",
    title = "Dominance strength in beluga samples by season"
  ) +
  theme_classic() +
  theme(legend.position = "none")

beluga_dominance_plot_season

ringed_dominance_plot_season <- ggplot(
  ringed_dominant_asv_per_sample %>%
    dplyr::filter(
      !is.na(season),
      trimws(as.character(season)) != ""
    ),
  aes(x = season, y = dominance, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    colour = "black"
  ) +
  geom_jitter(
    width = 0.15,
    shape = 21,
    colour = "black",
    size = 2,
    alpha = 0.8,
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Season",
    y = "Reads represented by dominant ASV",
    title = "Dominance strength in ringed seal samples by season"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ringed_dominance_plot_season

bearded_dominance_plot_season <- ggplot(
  bearded_dominant_asv_per_sample %>%
    dplyr::filter(
      !is.na(season),
      trimws(as.character(season)) != ""
    ),
  aes(x = season, y = dominance, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    colour = "black"
  ) +
  geom_jitter(
    width = 0.15,
    shape = 21,
    colour = "black",
    size = 2,
    alpha = 0.8,
    stroke = 0.35
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Season",
    y = "Reads represented by dominant ASV",
    title = "Dominance strength in bearded seal samples by season"
  ) +
  theme_classic() +
  theme(legend.position = "none")

bearded_dominance_plot_season

dominance_strength_plots_season <- (
  bearded_dominance_plot_season+
  beluga_dominance_plot_season +
    ringed_dominance_plot_season
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15)
  )

dominance_strength_plots_season

ggsave(
  "Deliverables/ubiome/dominant_asv_strength_by-season_each-host.png",
  plot = dominance_strength_plots_season,
  width = 18,
  height = 6,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# END DOMINANT-ASV SECTION
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Phylogenetics
# ------------------------------------------------------------------------------

# Faith

### biological interpretation and conclusion

# ------------------------------------------------------------------------------
# Information
# ------------------------------------------------------------------------------
# BY PREDATOR
# ------------------------------------------------------------------------------

# Kruskal-Wallis test
shan_test_host <- kruskal.test(
  Shannon ~ Predator,
  data = metadata
)

shan_test_host

# Dunn test to determine which host pairs differ
shan_dunn_results_host <- metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(Predator)
  ) %>%
  rstatix::dunn_test(
    Shannon ~ Predator,
    p.adjust.method = "hochberg"
  )

shan_dunn_results_host

# Add x- and y-coordinates for comparison brackets
shan_dunn_plot_host <- shan_dunn_results_host %>%
  rstatix::add_xy_position(x = "Predator")

shan_dunn_plot_host

# Box-and-whisker plot by host, with individual samples overlaid
shan_plot_host <- ggplot(
  metadata,
  aes(x = Predator, y = Shannon, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  ggpubr::stat_pvalue_manual(
    shan_dunn_plot_host,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Host",
    y = "Shannon entropy per sample",
    title = "Shannon entropy by host"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

shan_plot_host

# Gets global and by-host summary data
summary(metadata$Shannon)

metadata %>%
  dplyr::group_by(Predator) %>%
  dplyr::summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(median_shannon))

# ------------------------------------------------------------------------------
# Separate into individual hosts
# ------------------------------------------------------------------------------

# These match the naming convention used in your observed-richness section.
# You can retain them here, even if you already created these objects above.

beluga_metadata <- metadata %>%
  dplyr::filter(Predator == "beluga whale")

ringed_metadata <- metadata %>%
  dplyr::filter(Predator == "ringed seal")

bearded_metadata <- metadata %>%
  dplyr::filter(Predator == "bearded seal")

# Format age metadata for the seal-only age analysis
metadata_age <- metadata_age %>%
  dplyr::mutate(
    Age_Group = dplyr::case_when(
      is.na(Age_Group) |
        trimws(as.character(Age_Group)) %in% c("", "Pending") ~ "Unknown",
      TRUE ~ as.character(Age_Group)
    )
  )

ringed_metadata_age <- metadata_age %>%
  dplyr::filter(Predator == "ringed seal")

bearded_metadata_age <- metadata_age %>%
  dplyr::filter(Predator == "bearded seal")

# ------------------------------------------------------------------------------
# BY SEASON
# ------------------------------------------------------------------------------

# Exclude NA, blank, or whitespace-only seasonal values separately for each host.

beluga_metadata_season <- beluga_metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    season = droplevels(factor(season))
  )

ringed_metadata_season <- ringed_metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    season = droplevels(factor(season))
  )

bearded_metadata_season <- bearded_metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    season = droplevels(factor(season))
  )

# Check seasonal sample sizes
beluga_metadata_season %>%
  dplyr::count(season, name = "n")

ringed_metadata_season %>%
  dplyr::count(season, name = "n")

bearded_metadata_season %>%
  dplyr::count(season, name = "n")

# BELUGA: Shannon entropy by season
beluga_shan_test_season <- kruskal.test(
  Shannon ~ season,
  data = beluga_metadata_season
)

beluga_shan_test_season

beluga_shan_dunn_results_season <- beluga_metadata_season %>%
  rstatix::dunn_test(
    Shannon ~ season,
    p.adjust.method = "hochberg"
  )

beluga_shan_dunn_results_season

beluga_shan_dunn_plot_season <- beluga_shan_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

beluga_shan_dunn_plot_season

beluga_shan_plot_season <- ggplot(
  beluga_metadata_season,
  aes(x = season, y = Shannon, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  ggpubr::stat_pvalue_manual(
    beluga_shan_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Season",
    y = "Shannon entropy per beluga sample",
    title = "Beluga Shannon entropy by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

beluga_shan_plot_season

# RINGED SEAL: Shannon entropy by season
ringed_shan_test_season <- kruskal.test(
  Shannon ~ season,
  data = ringed_metadata_season
)

ringed_shan_test_season

ringed_shan_dunn_results_season <- ringed_metadata_season %>%
  rstatix::dunn_test(
    Shannon ~ season,
    p.adjust.method = "hochberg"
  )

ringed_shan_dunn_results_season

ringed_shan_dunn_plot_season <- ringed_shan_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

ringed_shan_dunn_plot_season

ringed_shan_plot_season <- ggplot(
  ringed_metadata_season,
  aes(x = season, y = Shannon, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  ggpubr::stat_pvalue_manual(
    ringed_shan_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Season",
    y = "Shannon entropy per ringed seal sample",
    title = "Ringed seal Shannon entropy by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ringed_shan_plot_season

# BEARDED SEAL: Shannon entropy by season
bearded_shan_test_season <- kruskal.test(
  Shannon ~ season,
  data = bearded_metadata_season
)

bearded_shan_test_season

bearded_shan_dunn_results_season <- bearded_metadata_season %>%
  rstatix::dunn_test(
    Shannon ~ season,
    p.adjust.method = "hochberg"
  )

bearded_shan_dunn_results_season

bearded_shan_dunn_plot_season <- bearded_shan_dunn_results_season %>%
  rstatix::add_xy_position(x = "season")

bearded_shan_dunn_plot_season

bearded_shan_plot_season <- ggplot(
  bearded_metadata_season,
  aes(x = season, y = Shannon, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  ggpubr::stat_pvalue_manual(
    bearded_shan_dunn_plot_season,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Season",
    y = "Shannon entropy per bearded seal sample",
    title = "Bearded seal Shannon entropy by season"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

bearded_shan_plot_season

# Summary statistics by host and season
shan_summary_season <- dplyr::bind_rows(
  beluga_metadata_season,
  ringed_metadata_season,
  bearded_metadata_season
) %>%
  dplyr::group_by(Predator, season) %>%
  dplyr::summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Predator, dplyr::desc(median_shannon))

shan_summary_season

# ------------------------------------------------------------------------------
# BY LOCALE
# ------------------------------------------------------------------------------

# Exclude NA, blank, or whitespace-only Locale values separately for each host.

beluga_metadata_locale <- beluga_metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    Locale = droplevels(factor(Locale))
  )

ringed_metadata_locale <- ringed_metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    Locale = droplevels(factor(Locale))
  )

bearded_metadata_locale <- bearded_metadata %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    Locale = droplevels(factor(Locale))
  )

# Check locale sample sizes
beluga_metadata_locale %>%
  dplyr::count(Locale, name = "n")

ringed_metadata_locale %>%
  dplyr::count(Locale, name = "n")

bearded_metadata_locale %>%
  dplyr::count(Locale, name = "n")

# BELUGA: Shannon entropy by locale
beluga_shan_test_locale <- kruskal.test(
  Shannon ~ Locale,
  data = beluga_metadata_locale
)

beluga_shan_test_locale

beluga_shan_dunn_results_locale <- beluga_metadata_locale %>%
  rstatix::dunn_test(
    Shannon ~ Locale,
    p.adjust.method = "hochberg"
  )

beluga_shan_dunn_results_locale

beluga_shan_dunn_plot_locale <- beluga_shan_dunn_results_locale %>%
  rstatix::add_xy_position(x = "Locale")

beluga_shan_dunn_plot_locale

beluga_shan_plot_locale <- ggplot(
  beluga_metadata_locale,
  aes(x = Locale, y = Shannon, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  ggpubr::stat_pvalue_manual(
    beluga_shan_dunn_plot_locale,
    label = "p.adj.signif",
    hide.ns = TRUE,
    tip.length = 0.01
  ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Locale",
    y = "Shannon entropy per beluga sample",
    title = "Beluga Shannon entropy by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

beluga_shan_plot_locale

# RINGED SEAL: Shannon entropy by locale
# ringed_shan_test_locale <- kruskal.test(
#   Shannon ~ Locale,
#   data = ringed_metadata_locale
# )
# 
# ringed_shan_test_locale
# 
# ringed_shan_dunn_results_locale <- ringed_metadata_locale %>%
#   rstatix::dunn_test(
#     Shannon ~ Locale,
#     p.adjust.method = "hochberg"
#   )
# 
# ringed_shan_dunn_results_locale
# 
# ringed_shan_dunn_plot_locale <- ringed_shan_dunn_results_locale %>%
#   rstatix::add_xy_position(x = "Locale")
# 
# ringed_shan_dunn_plot_locale

ringed_shan_plot_locale <- ggplot(
  ringed_metadata_locale,
  aes(x = Locale, y = Shannon, fill = Predator)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8,
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
  # ggpubr::stat_pvalue_manual(
  #   ringed_shan_dunn_plot_locale,
  #   label = "p.adj.signif",
  #   hide.ns = TRUE,
  #   tip.length = 0.01
  # ) +
  scale_fill_manual(
    values = host_colors,
    drop = FALSE,
    na.value = "grey70"
  ) +
  labs(
    x = "Locale",
    y = "Shannon entropy per ringed seal sample",
    title = "Ringed seal Shannon entropy by locale"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_classic() +
  theme(legend.position = "none")

ringed_shan_plot_locale

# BEARDED SEAL: Shannon entropy by locale
#
# All bearded samples may be Arctic. A Kruskal-Wallis/Dunn comparison requires
# at least two locale groups, so check the number of represented locale levels.

n_bearded_locales <- dplyr::n_distinct(bearded_metadata_locale$Locale)

if (n_bearded_locales >= 2) {
  
  bearded_shan_test_locale <- kruskal.test(
    Shannon ~ Locale,
    data = bearded_metadata_locale
  )
  
  bearded_shan_test_locale
  
  bearded_shan_dunn_results_locale <- bearded_metadata_locale %>%
    rstatix::dunn_test(
      Shannon ~ Locale,
      p.adjust.method = "hochberg"
    )
  
  bearded_shan_dunn_results_locale
  
  bearded_shan_dunn_plot_locale <- bearded_shan_dunn_results_locale %>%
    rstatix::add_xy_position(x = "Locale")
  
  bearded_shan_plot_locale <- ggplot(
    bearded_metadata_locale,
    aes(x = Locale, y = Shannon, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    ggpubr::stat_pvalue_manual(
      bearded_shan_dunn_plot_locale,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Locale",
      y = "Shannon entropy per bearded seal sample",
      title = "Bearded seal Shannon entropy by locale"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  bearded_shan_plot_locale
  
} else {
  
  message(
    "Bearded seal locale test skipped: only one Locale is represented ",
    "after dropping missing values."
  )
  
  bearded_shan_plot_locale <- ggplot(
    bearded_metadata_locale,
    aes(x = Locale, y = Shannon, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Locale",
      y = "Shannon entropy per bearded seal sample",
      title = "Bearded seal Shannon entropy by locale"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

bearded_shan_plot_locale

# Summary statistics by host and locale
shan_summary_locale <- dplyr::bind_rows(
  beluga_metadata_locale,
  ringed_metadata_locale,
  bearded_metadata_locale
) %>%
  dplyr::group_by(Predator, Locale) %>%
  dplyr::summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Predator, dplyr::desc(median_shannon))

shan_summary_locale

# ------------------------------------------------------------------------------
# BY AGE GROUP
# ------------------------------------------------------------------------------

# Age analyses include ringed and bearded seals only.
# Drop NA / blank / Pending / Unknown ages. Unknown is missing information,
# not a biological age group.

ringed_metadata_age_known <- ringed_metadata_age %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(Age_Group),
    trimws(as.character(Age_Group)) != "",
    Age_Group %in% c("Non-pup", "Pup")
  ) %>%
  dplyr::mutate(
    Age_Group = droplevels(
      factor(Age_Group, levels = c("Non-pup", "Pup"))
    )
  )

bearded_metadata_age_known <- bearded_metadata_age %>%
  dplyr::filter(
    !is.na(Shannon),
    !is.na(Age_Group),
    trimws(as.character(Age_Group)) != "",
    Age_Group %in% c("Non-pup", "Pup")
  ) %>%
  dplyr::mutate(
    Age_Group = droplevels(
      factor(Age_Group, levels = c("Non-pup", "Pup"))
    )
  )

# Check known-age group sample sizes
ringed_metadata_age_known %>%
  dplyr::count(Age_Group, name = "n")

bearded_metadata_age_known %>%
  dplyr::count(Age_Group, name = "n")

# RINGED SEAL: Shannon entropy by age group
n_ringed_age_groups <- dplyr::n_distinct(
  ringed_metadata_age_known$Age_Group
)
  
  ringed_shan_test_age <- kruskal.test(
    Shannon ~ Age_Group,
    data = ringed_metadata_age_known
  )
  
  ringed_shan_test_age
  
  ringed_shan_dunn_results_age <- ringed_metadata_age_known %>%
    rstatix::dunn_test(
      Shannon ~ Age_Group,
      p.adjust.method = "hochberg"
    )
  
  ringed_shan_dunn_results_age
  
  ringed_shan_dunn_plot_age <- ringed_shan_dunn_results_age %>%
    rstatix::add_xy_position(x = "Age_Group")
  
  ringed_shan_plot_age <- ggplot(
    ringed_metadata_age_known,
    aes(x = Age_Group, y = Shannon, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    ggpubr::stat_pvalue_manual(
      ringed_shan_dunn_plot_age,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Age group",
      y = "Shannon entropy per ringed seal sample",
      title = "Ringed seal Shannon entropy by age group"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  ringed_shan_plot_age


  ringed_shan_plot_age <- ggplot(
    ringed_metadata_age_known,
    aes(x = Age_Group, y = Shannon, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Age group",
      y = "Shannon entropy per ringed seal sample",
      title = "Ringed seal Shannon entropy by age group"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")

ringed_shan_plot_age

# BEARDED SEAL: Shannon entropy by age group
n_bearded_age_groups <- dplyr::n_distinct(
  bearded_metadata_age_known$Age_Group
)

  bearded_shan_test_age <- kruskal.test(
    Shannon ~ Age_Group,
    data = bearded_metadata_age_known
  )
  
  bearded_shan_test_age
  
  bearded_shan_dunn_results_age <- bearded_metadata_age_known %>%
    rstatix::dunn_test(
      Shannon ~ Age_Group,
      p.adjust.method = "hochberg"
    )
  
  bearded_shan_dunn_results_age
  
  bearded_shan_dunn_plot_age <- bearded_shan_dunn_results_age %>%
    rstatix::add_xy_position(x = "Age_Group")
  
  bearded_shan_plot_age <- ggplot(
    bearded_metadata_age_known,
    aes(x = Age_Group, y = Shannon, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    ggpubr::stat_pvalue_manual(
      bearded_shan_dunn_plot_age,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Age group",
      y = "Shannon entropy per bearded seal sample",
      title = "Bearded seal Shannon entropy by age group"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  bearded_shan_plot_age
  
  
  bearded_shan_plot_age <- ggplot(
    bearded_metadata_age_known,
    aes(x = Age_Group, y = Shannon, fill = Predator)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      size = 2,
      alpha = 0.8,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = host_colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = "Age group",
      y = "Shannon entropy per bearded seal sample",
      title = "Bearded seal Shannon entropy by age group"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_classic() +
    theme(legend.position = "none")

bearded_shan_plot_age

# Summary statistics by host and known age group
shan_summary_age <- dplyr::bind_rows(
  ringed_metadata_age_known,
  bearded_metadata_age_known
) %>%
  dplyr::group_by(Predator, Age_Group) %>%
  dplyr::summarise(
    n = sum(!is.na(Shannon)),
    median_shannon = median(Shannon, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    IQR_shannon = IQR(Shannon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Predator, dplyr::desc(median_shannon))

shan_summary_age

# ------------------------------------------------------------------------------
# SHANNON PLOTS
# ------------------------------------------------------------------------------

# Host-level Shannon plot
shan_plot_host

# Seasonal Shannon plots: one for each host
shan_plots_season <- (
  bearded_shan_plot_season+
  beluga_shan_plot_season +
    ringed_shan_plot_season
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

shan_plots_season

# ggsave(
#   "Deliverables/ubiome/alphadiv/SHAN_ubiome-majorclass_by-season.png",
#   plot = shan_plots_season,
#   width = 18,
#   height = 6,
#   units = "in",
#   dpi = 300
# )

# Locale Shannon plots: one for each host
shan_plots_locale <- (
  bearded_shan_plot_locale+
  beluga_shan_plot_locale +
    ringed_shan_plot_locale
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

shan_plots_locale

# ggsave(
#   "Deliverables/ubiome/alphadiv/SHAN_ubiome-majorclass_by-locale.png",
#   plot = shan_plots_locale,
#   width = 18,
#   height = 6,
#   units = "in",
#   dpi = 300
# )

# Age-group Shannon plots: ringed and bearded seals only
shan_plots_age <- (
  bearded_shan_plot_age+
  ringed_shan_plot_age
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

shan_plots_age

# ggsave(
#   "Deliverables/ubiome/alphadiv/SHAN_ubiome-majorclass_by-age.png",
#   plot = shan_plots_age,
#   width = 12,
#   height = 6,
#   units = "in",
#   dpi = 300
# )

# Combined Shannon figure
shan_plots <- (
  shan_plot_host /
    shan_plots_season /
    shan_plots_locale /
    shan_plots_age
) &
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

shan_plots

ggsave(
  "Deliverables/ubiome/alphadiv/SHAN_ubiome-majorclass.png",
  plot = shan_plots,
  width = 20,
  height = 26,
  units = "in",
  dpi = 300
)

# interpretation of all metrics collectively

