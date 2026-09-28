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
load("./Scripts/12s/rdata/ps.12s.dietdat.Rdata")

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

# leftjoin diet alphadiv info to metadata

metadata_12s <- metadata_12s%>%
  dplyr::rename(diet_observed = Observed)%>%
  dplyr::rename(diet_bp = Berger_Parker)%>%
  dplyr::rename(diet_shan = Shannon)

metadata <- metadata %>%
  dplyr::left_join(
    metadata_12s %>%
      dplyr::select(LabID, diet_observed, diet_bp, diet_shan),
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
norm # p-value = 0.013 (significant - not normally distributed)

# test homogeneity of variance FOR HOST
#?leveneTest()
var <- leveneTest(Observed ~ Predator, metadata)
var # p-value = 0.02118 (significant - non-homonegeous variance)

var <- leveneTest(Observed ~ Month, metadata)
var # p-value = 0.3756 (non-significant - variance is homogeneous)

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

# p-value = 0.8556; though there is a slightly negative trend, there is no evidence
# that the read depth is highly correlated with the observed richness
### interpretation = it's okay not to rarefy (which I don't want to do anyway!)



# BY PREDATOR----------------------------------------------------------------------------
# Kruskall-Wallis test
obs_test_host <- kruskal.test(Observed ~ Predator, data = metadata)
obs_test_host #  p-value = 6.034e-05 (significant)

# Dunn test to determine which groups are significantly different
#?dunn_test()
obs_dunn_results_host <- metadata %>%
  dunn_test(
    Observed ~ Predator,
    # p.adjust.method = "hochberg"
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

# DIET TESTS-------------------------------------------------------------------------------------

# function to check for correlations
run_spearman <- function(data, host, diet_metric, label) {
  
  dat <- data %>%
    dplyr::filter(
      Predator == host,
      !is.na(Observed),
      !is.na(.data[[diet_metric]])
    )
  
  test <- cor.test(
    dat$Observed,
    dat[[diet_metric]],
    method = "spearman",
    exact = FALSE
  )
  
  tibble::tibble(
    Predator = host,
    diet_metric = label,
    n = nrow(dat),
    rho = unname(test$estimate),
    p_raw = test$p.value
  )
}

# spearman tests
observed_diet_tests <- dplyr::bind_rows(
  
  run_spearman(
    metadata,
    host = "beluga whale",
    diet_metric = "diet_observed",
    label = "Diet observed richness"
  ),
  run_spearman(
    metadata,
    host = "bearded seal",
    diet_metric = "diet_observed",
    label = "Diet observed richness"
  ),
  run_spearman(
    metadata,
    host = "ringed seal",
    diet_metric = "diet_observed",
    label = "Diet observed richness"
  ),
  
  run_spearman(
    metadata,
    host = "beluga whale",
    diet_metric = "diet_bp",
    label = "Diet Berger-Parker dominance"
  ),
  run_spearman(
    metadata,
    host = "bearded seal",
    diet_metric = "diet_bp",
    label = "Diet Berger-Parker dominance"
  ),
  run_spearman(
    metadata,
    host = "ringed seal",
    diet_metric = "diet_bp",
    label = "Diet Berger-Parker dominance"
  ),
  
  run_spearman(
    metadata,
    host = "beluga whale",
    diet_metric = "diet_shan",
    label = "Diet Shannon diversity"
  ),
  run_spearman(
    metadata,
    host = "bearded seal",
    diet_metric = "diet_shan",
    label = "Diet Shannon diversity"
  ),
  run_spearman(
    metadata,
    host = "ringed seal",
    diet_metric = "diet_shan",
    label = "Diet Shannon diversity"
  )
)

observed_diet_tests

## beluga diet observed richness, Berger-Parker dominance, and Shannon diversity are significant 
### this is based on the raw p value -> a correction may be smart here, as I tested many correlations.


# tests appropriateness of GAM for geom_smooth on plots for richness
metadata %>%
  dplyr::filter(
    Predator %in% names(host_colors),
    !is.na(Observed),
    !is.na(diet_observed)
  ) %>%
  dplyr::group_by(Predator) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_unique_diet_richness = dplyr::n_distinct(diet_observed),
    min_diet_richness = min(diet_observed),
    max_diet_richness = max(diet_observed),
    .groups = "drop"
  )
## too few unique x values for a GAM. Try GLM or GLMM instead.

# richness plots
plot_observed_vs_diet_richness <- metadata %>%
  dplyr::filter(
    Predator %in% names(host_colors),
    !is.na(Observed),
    !is.na(diet_observed)
  ) %>%
  ggplot(
    aes(
      x = diet_observed,
      y = Observed,
      color = Predator
    )
  ) +
  geom_point(
    size = 2.3,
    alpha = 0.8
  ) +
  facet_wrap(~ Predator, nrow = 1) +
  scale_color_manual(values = host_colors) +
  theme_classic() +
  labs(
    x = "Diet observed ASV richness (12S)",
    y = "Microbiome observed ASV richness",
    title = "Microbiome richness versus diet richness"
  ) +
  theme(legend.position = "none")

plot_observed_vs_diet_richness

# dominance plots
plot_observed_vs_diet_bp <- metadata %>%
  dplyr::filter(
    Predator %in% names(host_colors),
    !is.na(Observed),
    !is.na(diet_bp)
  ) %>%
  ggplot(
    aes(
      x = diet_bp,
      y = Observed,
      color = Predator
    )
  ) +
  geom_point(
    size = 2.3,
    alpha = 0.8
  ) +
  facet_wrap(~ Predator, nrow = 1) +
  scale_color_manual(values = host_colors) +
  theme_classic() +
  labs(
    x = "Diet Berger–Parker dominance (12S)",
    y = "Microbiome observed ASV richness",
    title = "Microbiome richness versus diet dominance"
  ) +
  theme(legend.position = "none")

plot_observed_vs_diet_bp

# shannon plots
plot_observed_vs_diet_shannon <- metadata %>%
  dplyr::filter(
    Predator %in% names(host_colors),
    !is.na(Observed),
    !is.na(diet_shan)
  ) %>%
  ggplot(
    aes(
      x = diet_shan,
      y = Observed,
      color = Predator
    )
  ) +
  geom_point(
    size = 2.3,
    alpha = 0.8
  ) +
  facet_wrap(~ Predator, nrow = 1) +
  scale_color_manual(values = host_colors) +
  theme_classic() +
  labs(
    x = "Diet Shannon diversity (12S)",
    y = "Microbiome observed ASV richness",
    title = "Microbiome richness versus diet Shannon diversity"
  ) +
  theme(legend.position = "none")

plot_observed_vs_diet_shannon

