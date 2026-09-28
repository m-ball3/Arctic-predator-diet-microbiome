# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS ZETA DIVERSITY METRICS ON UBIOME DATA
# THIS IS THE SECOND STATISTICAL INVESTIGATION 
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
# ------------------------------------------------------------------------------

# install.packages(
#   "rbiom",
#   type = "binary",
#   repos = "https://cran.rstudio.com"
# )

# remotes::install_version(
#   package = "zetadiv",
#   version = "1.2.1",
#   repos = "https://cran.r-project.org"
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
library(zetadiv)

??zetadiv()
getwd()
load("./Scripts/12s/rdata/ps.12s.dietdat.Rdata")
load("./Scripts/ubiome/rdata/ubiome.ps.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.norep), "matrix")
head(rownames(asv_mat))

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
  dplyr::filter(!row_id %in% c("WADE-003-150-A", "WADE-003-150-B",
                               "neg1","mock1-782D", "mock2-783D", "neg2","neg2-2" )) %>% 
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


# Samples in metadata but not in asv mat
setdiff(rownames(metadata), rownames(asv_mat))

# Samples in asv mat but not in metadata
setdiff(rownames(asv_mat), rownames(metadata))

# formatting metadata for age class (only bearded and ringed seals)
metadata_age <- metadata %>%
  dplyr::filter(Predator %in% c("ringed seal", "bearded seal")) %>%
  dplyr::mutate(
    Age_Group = dplyr::case_when(
      Age_Group %in% c("", "Pending") ~ "Unknown",
      TRUE ~ Age_Group
    )
  )
metadata_age <- metadata_age %>%
  dplyr::filter(!Age_Group == "Unknown")

# formatting asv_mat for age class (only bearded and ringed seals)
asv_mat_age <- asv_mat %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

# ------------------------------------------------------------------------------
# Calculates zeta div; creates table
# ------------------------------------------------------------------------------

# creates a presence - absence table
## 1 = the ASV was present in the sample
## 0 = the ASV was absence in the sample
asv_pa <- (asv_mat > 0) * 1
summary(asv_pa) #sanity check

# FOR HOST
# ensures correct ordering of the metadata table
sample_ids_by_predator <- split(
  rownames(metadata),
  metadata$Predator
)

# creates a list object with the presence - absence infromation, 
# the predator assignment, and the number of samples per predator
zeta.pred <- purrr::imap(sample_ids_by_predator, function(sample_ids, predator_name) {
  
  pa_group <- as.data.frame(asv_pa[sample_ids, , drop = FALSE])
  n_sites <- nrow(pa_group)
  
  out <- zetadiv::Zeta.decline.ex( # calculates expected numbers of taxa shared among combinations of samples at successive zeta orders
    data.spec = pa_group,
    orders = seq_len(n_sites),
    plot = FALSE
  )
  
  list(
    predator = predator_name,
    n_samples = n_sites,
    result = out
  )
})

names(zeta.pred) # sanity check

# creates table of number of samples per predator, zeta order and number, and retention
zeta_summary_host <- imap_dfr(zeta.pred, function(x, predator_name) {
  
  zeta_vals <- x$result$zeta.val
  
  tibble(
    Predator = predator_name,
    n_samples = x$n_samples,
    order = seq_along(zeta_vals),
    zeta = zeta_vals,
    retention = c(NA, zeta_vals[-1] / zeta_vals[-length(zeta_vals)])
  )
})

zeta_summary_host

# the fraction of ASVs retained as one more sample is included. A steep decline 
# indicates substantial among-sample turnover; sustained retention at higher 
# orders indicates a relatively consistent core assemblage within that predator group.


# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

# Zeta order (i) = number of samples / assemblages included in a comparison.
# It ranges from 1 to the total number of samples.

# Zeta_i = mean number of taxa shared among all combinations of i samples.
# Example: zeta_10 is the mean number of taxa shared among each
# combination of 10 samples—not necessarily only samples 1–10.

# Zeta_(i - 1) = mean number of taxa shared among all combinations of
# i - 1 samples.
# Example: zeta_(10-1) is the mean number of taxa shared among each
# combination of 9 samples—not specifically samples 1–9.

# shows how the number of ASVs shared across multiple samples 
# declines as you require an ASV to occur in more samples. 
p1_host <- ggplot(
  zeta_summary_host,
  aes(x = order, y = zeta, colour = Predator, group = Predator)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  labs(
    x = "",
    y = "Mean shared ASVs",
    colour = "Predator",
    title = "Mean Shared ASVs by Zeta Order"
  ) +
  theme_classic()

p1_host # interp: 
  

# y = z component ratio = Zeta_i / Zeta_i-1
## the probability that a species shared by i-1 sites is also 
## found to be shared by i sites can be expressed as the z component ratio 

## zeta order - the number of sites (aka the number of samples in each group)

# shows the rate at which shared ASVs are retained as another sample is added.
p2_host <- ggplot(
  filter(zeta_summary_host, order > 1),
  aes(x = order, y = retention, colour = Predator, group = Predator)
) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    x = "",
    y = expression(zeta[i] / zeta[i-1]),
    colour = "Predator",
    title = "Component Ratio by Zeta Order"
  ) +
  theme_classic()

p2_host # interp
# At each order, it is the proportion of ASVs shared across
# i−1 samples that remain shared when requiring presence in one additional sample.

# A value near 1: nearly all currently shared ASVs persist into the next order; low additional turnover.

# A value near 0: virtually none remain when another sample is added; high turnover.

# The dashed line at 1 is the theoretical maximum.

# ==============================================================================
# ZETA DIVERSITY WITHIN EACH HOST
#
# Analyses:
#   - Season within beluga whale, ringed seal, and bearded seal samples
#   - Locale within beluga whale and ringed seal samples
#   - Age group within ringed and bearded seal samples
#
# Location is intentionally excluded.
#
# This code assumes:
#   - asv_pa has samples in rows and ASVs in columns.
#   - rownames(asv_pa) are LabID values.
#   - metadata contains LabID, Predator, season, and Locale.
#   - metadata_age contains LabID, Predator, and Age_Group.
# ==============================================================================

# ------------------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------------------

# Maximum zeta order shown/calculated in each subgroup.
# A common order is important for comparing lines among covariate groups.
# The code automatically uses smaller maximum orders for groups with fewer samples.
max_zeta_order <- 15

# Number of random site combinations sampled per zeta order for the
# Monte-Carlo zeta-decline calculation.
zeta_samples <- 1000

# ------------------------------------------------------------------------------ 
# FUNCTION: COMPUTE ZETA DECLINE FOR EACH LEVEL OF A COVARIATE WITHIN ONE HOST
# ------------------------------------------------------------------------------

run_zeta_within_host <- function(
    metadata_df,
    pa_matrix,
    host_name,
    covariate,
    sample_id = "LabID",
    host_var = "Predator",
    max_order = 15,
    sam = 1000,
    exclude_values = character(),
    min_n_per_group = 2
) {
  
  # Convert inputs to stable types.
  metadata_df <- metadata_df %>%
    dplyr::mutate(
      "{sample_id}" := as.character(.data[[sample_id]]),
      "{host_var}" := as.character(.data[[host_var]])
    ) %>%
    dplyr::filter(
      .data[[host_var]] == host_name,
      .data[[sample_id]] %in% rownames(pa_matrix),
      !is.na(.data[[covariate]]),
      trimws(as.character(.data[[covariate]])) != "",
      !as.character(.data[[covariate]]) %in% exclude_values
    )
  
  # Preserve presence-absence matrix order, then arrange metadata to match.
  pa_sub_all <- pa_matrix[
    rownames(pa_matrix) %in% metadata_df[[sample_id]],
    ,
    drop = FALSE
  ]
  
  metadata_df <- metadata_df %>%
    dplyr::slice(
      match(
        rownames(pa_sub_all),
        .data[[sample_id]]
      )
    ) %>%
    dplyr::mutate(
      "{covariate}" := droplevels(factor(.data[[covariate]]))
    )
  
  stopifnot(
    identical(
      as.character(metadata_df[[sample_id]]),
      rownames(pa_sub_all)
    )
  )
  
  # Split the matching sample IDs by covariate level.
  sample_ids_by_group <- split(
    metadata_df[[sample_id]],
    metadata_df[[covariate]]
  )
  
  # Calculate zeta decline separately for every level.
  zeta_results <- purrr::imap(
    sample_ids_by_group,
    function(sample_ids, group_name) {
      
      pa_group <- as.data.frame(
        pa_sub_all[sample_ids, , drop = FALSE]
      )
      
      n_sites <- nrow(pa_group)
      
      # Skip groups too small to calculate a decline.
      if (n_sites < min_n_per_group) {
        return(
          list(
            group = group_name,
            n_samples = n_sites,
            result = NULL,
            status = "Skipped: fewer than 2 samples"
          )
        )
      }
      
      # Cannot calculate a zeta order above the group sample count.
      orders_to_use <- seq_len(min(max_order, n_sites))
      
      set.seed(123)
      
      # Monte Carlo calculation.
      zeta_out <- zetadiv::Zeta.decline.mc(
        data.spec = pa_group,
        orders = orders_to_use,
        sam = sam,
        plot = FALSE,
        silent = TRUE
      )
      
      list(
        group = group_name,
        n_samples = n_sites,
        result = zeta_out,
        status = "OK"
      )
    }
  )
  
  # Convert results into one tidy output data frame.
  zeta_summary <- purrr::imap_dfr(
    zeta_results,
    function(x, group_name) {
      
      if (is.null(x$result)) {
        return(
          tibble::tibble(
            Predator = host_name,
            covariate = covariate,
            group = group_name,
            n_samples = x$n_samples,
            order = NA_integer_,
            zeta = NA_real_,
            retention = NA_real_,
            status = x$status
          )
        )
      }
      
      zeta_vals <- x$result$zeta.val
      
      tibble::tibble(
        Predator = host_name,
        covariate = covariate,
        group = group_name,
        n_samples = x$n_samples,
        order = seq_along(zeta_vals),
        zeta = zeta_vals,
        retention = c(
          NA_real_,
          zeta_vals[-1] / zeta_vals[-length(zeta_vals)]
        ),
        status = x$status
      )
    }
  )
  
  zeta_summary
}

# ------------------------------------------------------------------------------
# SEASON WITHIN EACH HOST
# ------------------------------------------------------------------------------

# Remove NA, empty, and whitespace-only season records for each host separately.

zeta_beluga_season <- run_zeta_within_host(
  metadata_df = metadata,
  pa_matrix = asv_pa,
  host_name = "beluga whale",
  covariate = "season",
  max_order = max_zeta_order,
  sam = zeta_samples
)

zeta_ringed_season <- run_zeta_within_host(
  metadata_df = metadata,
  pa_matrix = asv_pa,
  host_name = "ringed seal",
  covariate = "season",
  max_order = max_zeta_order,
  sam = zeta_samples
)

zeta_bearded_season <- run_zeta_within_host(
  metadata_df = metadata,
  pa_matrix = asv_pa,
  host_name = "bearded seal",
  covariate = "season",
  max_order = max_zeta_order,
  sam = zeta_samples
)

zeta_summary_season <- dplyr::bind_rows(
  zeta_beluga_season,
  zeta_ringed_season,
  zeta_bearded_season
)

zeta_summary_season

# Check sample sizes per season within each host.
zeta_summary_season %>%
  dplyr::distinct(Predator, group, n_samples, status) %>%
  dplyr::arrange(Predator, group)

# ------------------------------------------------------------------------------
# LOCALE WITHIN EACH HOST
# ------------------------------------------------------------------------------

# Bearded seals are not included because all bearded-seal samples occur in the
# Arctic, so there is no within-host locale comparison.

zeta_beluga_locale <- run_zeta_within_host(
  metadata_df = metadata,
  pa_matrix = asv_pa,
  host_name = "beluga whale",
  covariate = "Locale",
  max_order = max_zeta_order,
  sam = zeta_samples
)

zeta_ringed_locale <- run_zeta_within_host(
  metadata_df = metadata,
  pa_matrix = asv_pa,
  host_name = "ringed seal",
  covariate = "Locale",
  max_order = max_zeta_order,
  sam = zeta_samples
)

zeta_bearded_locale <- run_zeta_within_host(
  metadata_df = metadata,
  pa_matrix = asv_pa,
  host_name = "bearded seal",
  covariate = "Locale",
  max_order = max_zeta_order,
  sam = zeta_samples
)

zeta_summary_locale <- dplyr::bind_rows(
  zeta_bearded_locale,
  zeta_beluga_locale,
  zeta_ringed_locale
)

zeta_summary_locale

# Check sample sizes per locale within each host.
zeta_summary_locale %>%
  dplyr::distinct(Predator, group, n_samples, status) %>%
  dplyr::arrange(Predator, group)

# ------------------------------------------------------------------------------
# AGE GROUP WITHIN EACH SEAL HOST
# ------------------------------------------------------------------------------

# Exclude missing-data categories from primary biological age-group comparisons.
# metadata_age should already contain only ringed and bearded seal samples.

zeta_ringed_age <- run_zeta_within_host(
  metadata_df = metadata_age,
  pa_matrix = asv_pa,
  host_name = "ringed seal",
  covariate = "Age_Group",
  max_order = max_zeta_order,
  sam = zeta_samples,
  exclude_values = c("", "Pending", "Unknown")
)

zeta_bearded_age <- run_zeta_within_host(
  metadata_df = metadata_age,
  pa_matrix = asv_pa,
  host_name = "bearded seal",
  covariate = "Age_Group",
  max_order = max_zeta_order,
  sam = zeta_samples,
  exclude_values = c("", "Pending", "Unknown")
)

zeta_summary_age <- dplyr::bind_rows(
  zeta_ringed_age,
  zeta_bearded_age
)

zeta_summary_age

# Check sample sizes per known age class within each seal host.
zeta_summary_age %>%
  dplyr::distinct(Predator, group, n_samples, status) %>%
  dplyr::arrange(Predator, group)

# ------------------------------------------------------------------------------
# FUNCTION: ZETA-DECLINE PLOT
# ------------------------------------------------------------------------------

plot_zeta_decline <- function(
    zeta_summary,
    title,
    x_label = "Zeta order (number of samples)",
    y_label = "Mean shared ASVs"
) {
  
  ggplot(
    zeta_summary %>%
      dplyr::filter(
        status == "OK",
        !is.na(order),
        !is.na(zeta)
      ),
    aes(
      x = order,
      y = zeta,
      colour = group,
      group = group
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = scales::breaks_pretty()
    ) +
    labs(
      x = x_label,
      y = y_label,
      colour = NULL,
      title = title
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------------------------
# FUNCTION: ZETA-RETENTION PLOT
# ------------------------------------------------------------------------------

plot_zeta_retention <- function(
    zeta_summary,
    title,
    x_label = "Zeta order (number of samples)"
) {
  
  ggplot(
    zeta_summary %>%
      dplyr::filter(
        status == "OK",
        order > 1,
        !is.na(retention)
      ),
    aes(
      x = order,
      y = retention,
      colour = group,
      group = group
    )
  ) +
    geom_hline(
      yintercept = 1,
      linetype = 2,
      colour = "grey60"
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = scales::breaks_pretty()
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = x_label,
      y = expression(zeta[i] / zeta[i-1]),
      colour = NULL,
      title = title
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------------------------
# SEASON PLOTS BY HOST
# ------------------------------------------------------------------------------

p1_beluga_season <- plot_zeta_decline(
  zeta_beluga_season,
  title = "Beluga whale zeta decline by season"
)

p2_beluga_season <- plot_zeta_retention(
  zeta_beluga_season,
  title = "Beluga whale ASV retention by season"
)

p1_ringed_season <- plot_zeta_decline(
  zeta_ringed_season,
  title = "Ringed seal zeta decline by season"
)

p2_ringed_season <- plot_zeta_retention(
  zeta_ringed_season,
  title = "Ringed seal ASV retention by season"
)

p1_bearded_season <- plot_zeta_decline(
  zeta_bearded_season,
  title = "Bearded seal zeta decline by season"
)

p2_bearded_season <- plot_zeta_retention(
  zeta_bearded_season,
  title = "Bearded seal ASV retention by season"
)

zeta_plots_season <- (
  (p1_beluga_season + p2_beluga_season) /
    (p1_ringed_season + p2_ringed_season) /
    (p1_bearded_season + p2_bearded_season)
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

zeta_plots_season

ggsave(
  "Deliverables/ubiome/zetadiv/ZETA_ubiome-majorclass_by-season_each-host.png",
  plot = zeta_plots_season,
  width = 16,
  height = 18,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# LOCALE PLOTS BY HOST
# ------------------------------------------------------------------------------

p1_beluga_locale <- plot_zeta_decline(
  zeta_beluga_locale,
  title = "Beluga whale zeta decline by locale"
)

p2_beluga_locale <- plot_zeta_retention(
  zeta_beluga_locale,
  title = "Beluga whale ASV retention by locale"
)

p1_ringed_locale <- plot_zeta_decline(
  zeta_ringed_locale,
  title = "Ringed seal zeta decline by locale"
)

p2_ringed_locale <- plot_zeta_retention(
  zeta_ringed_locale,
  title = "Ringed seal ASV retention by locale"
)

zeta_plots_locale <- (
  (p1_beluga_locale + p2_beluga_locale) /
    (p1_ringed_locale + p2_ringed_locale)
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

zeta_plots_locale

ggsave(
  "Deliverables/ubiome/zetadiv/ZETA_ubiome-majorclass_by-locale_each-host.png",
  plot = zeta_plots_locale,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# AGE-GROUP PLOTS BY HOST
# ------------------------------------------------------------------------------

p1_ringed_age <- plot_zeta_decline(
  zeta_ringed_age,
  title = "Ringed seal zeta decline by age group"
)

p2_ringed_age <- plot_zeta_retention(
  zeta_ringed_age,
  title = "Ringed seal ASV retention by age group"
)

p1_bearded_age <- plot_zeta_decline(
  zeta_bearded_age,
  title = "Bearded seal zeta decline by age group"
)

p2_bearded_age <- plot_zeta_retention(
  zeta_bearded_age,
  title = "Bearded seal ASV retention by age group"
)

zeta_plots_age <- (
  (p1_ringed_age + p2_ringed_age) /
    (p1_bearded_age + p2_bearded_age)
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

zeta_plots_age

ggsave(
  "Deliverables/ubiome/zetadiv/ZETA_ubiome-majorclass_by-age_each-seal-host.png",
  plot = zeta_plots_age,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300
)

# ------------------------------------------------------------------------------
# FINAL COMBINED ZETA FIGURE
# ------------------------------------------------------------------------------

all_zeta_plots <- (
  zeta_plots_season /
    zeta_plots_locale /
    zeta_plots_age
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

all_zeta_plots

ggsave(
  "Deliverables/ubiome/zetadiv/ZETA_ubiome-majorclass_by-host.png",
  plot = all_zeta_plots,
  width = 18,
  height = 40,
  units = "in",
  dpi = 300
)

# ==============================================================================
# BEGIN: CORE TAXA SECTION
# this code seeks to find a "core" microbiome across combinations of these samples
## a "core" microbiome, here is defined when: 
### An ASV was classified as a candidate host-specific core member if detected 
### in at least 80% of samples within that host and with mean relative abundance 
### at least 0.1% across those samples.
# ==============================================================================

# are there any samples in 100% of the samples tested?
colSums(asv_pa > 0) == nrow(asv_pa) # no ASVs present in 100% of the samples in asv_pa.

# creates a "core" ASV table, including:
#-------------------------------------------------------------------------------------------
# number of samples with the ASV,
# 
# occupancy proportion,
# 
# mean relative abundance across all samples,
# 
# mean relative abundance among samples in which it was detected,
# 
# total reads,
# 
# rank by occupancy and abundance,
# 
# core membership under your selected thresholds.
#-------------------------------------------------------------------------------------------

# function
make_core_table <- function(
    count_matrix,
    tax_table_df,
    prevalence_threshold = 0.80,
    detection_threshold = 0.001
) {
  
  count_matrix <- as.matrix(count_matrix)
  
  # Presence/absence, using > 0 reads as detection.
  pa_matrix <- (count_matrix > 0) * 1
  
  # Relative abundance within each sample.
  rel_matrix <- sweep(
    count_matrix,
    1,
    rowSums(count_matrix),
    "/"
  )
  
  # Safeguard against any zero-library samples.
  rel_matrix[!is.finite(rel_matrix)] <- 0
  
  taxa_table <- tax_table_df %>%
    tibble::rownames_to_column("category")
  
  core_table <- tibble::tibble(
    category = colnames(count_matrix),
    n_samples = nrow(count_matrix),
    n_present = colSums(pa_matrix),
    prevalence = colMeans(pa_matrix),
    total_reads = colSums(count_matrix),
    mean_rel_abundance_all = colMeans(rel_matrix),
    median_rel_abundance_all = apply(rel_matrix, 2, median),
    mean_rel_abundance_present = vapply(
      seq_len(ncol(rel_matrix)),
      function(j) {
        x <- rel_matrix[, j]
        if (any(x > 0)) mean(x[x > 0]) else 0
      },
      numeric(1)
    ),
    max_rel_abundance = apply(rel_matrix, 2, max)
  ) %>%
    dplyr::mutate(
      core_prevalence = prevalence >= prevalence_threshold,
      core_abundance = mean_rel_abundance_all >= detection_threshold,
      core_candidate = core_prevalence & core_abundance
    ) %>%
    dplyr::left_join(
      taxa_table,
      by = "category"
    ) %>%
    dplyr::arrange(
      dplyr::desc(core_candidate),
      dplyr::desc(prevalence),
      dplyr::desc(mean_rel_abundance_all)
    )
  
  core_table
}

# uses function
tax_df_full <- as.data.frame(
  phyloseq::tax_table(ps.major)
)

core_all_hosts <- make_core_table(
  count_matrix = asv_mat,
  tax_table_df = tax_df_full,
  prevalence_threshold = 0.80,
  detection_threshold = 0.001
)

core_all_hosts

#strict universal core
strict <- core_all_hosts %>%
  dplyr::filter(prevalence == 1)

# 80% 
in_80_percent <- core_all_hosts %>%
  dplyr::filter(core_candidate)


# CORE SEPARATED BY HOST -----------------------------------------------------------------
# function
make_host_core_table <- function(
    host_name,
    metadata_df,
    count_matrix,
    tax_table_df,
    prevalence_threshold = 0.80,
    detection_threshold = 0.001
) {
  
  host_ids <- metadata_df %>%
    dplyr::filter(Predator == host_name) %>%
    dplyr::pull(LabID) %>%
    intersect(rownames(count_matrix))
  
  host_counts <- count_matrix[
    host_ids,
    ,
    drop = FALSE
  ]
  
  make_core_table(
    count_matrix = host_counts,
    tax_table_df = tax_table_df,
    prevalence_threshold = prevalence_threshold,
    detection_threshold = detection_threshold
  ) %>%
    dplyr::mutate(
      Predator = host_name,
      .before = 1
    )
}

# use function
core_beluga <- make_host_core_table(
  host_name = "beluga whale",
  metadata_df = metadata,
  count_matrix = asv_mat,
  tax_table_df = tax_df_full,
  prevalence_threshold = 0.80,
  detection_threshold = 0.001
)

core_bearded <- make_host_core_table(
  host_name = "bearded seal",
  metadata_df = metadata,
  count_matrix = asv_mat,
  tax_table_df = tax_df_full,
  prevalence_threshold = 0.80,
  detection_threshold = 0.001
)

core_ringed <- make_host_core_table(
  host_name = "ringed seal",
  metadata_df = metadata,
  count_matrix = asv_mat,
  tax_table_df = tax_df_full,
  prevalence_threshold = 0.80,
  detection_threshold = 0.001
)

core_by_host <- dplyr::bind_rows(
  core_beluga,
  core_bearded,
  core_ringed
)

core_by_host %>%
  dplyr::filter(core_candidate) %>%
  dplyr::select(
    Predator,
    category,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    #Species,
    n_samples,
    n_present,
    prevalence,
    mean_rel_abundance_all,
    mean_rel_abundance_present,
    max_rel_abundance,
    total_reads
  ) %>%
  dplyr::arrange(
    Predator,
    dplyr::desc(prevalence),
    dplyr::desc(mean_rel_abundance_all)
  )


# HOST BY ZETA----------------------------------------------------------------------------
core_from_zeta_groups <- purrr::imap_dfr(
  sample_ids_by_predator,
  function(sample_ids, predator_name) {
    
    make_core_table(
      count_matrix = asv_mat[sample_ids, , drop = FALSE],
      tax_table_df = tax_df_full,
      prevalence_threshold = 0.80,
      detection_threshold = 0.001
    ) %>%
      dplyr::mutate(
        Predator = predator_name,
        .before = 1
      )
  }
)

# HELPFUL SUMMARY TABLE
core_summary_by_host <- core_from_zeta_groups %>%
  dplyr::group_by(Predator) %>%
  dplyr::summarise(
    n_samples = dplyr::first(n_samples),
    universal_ASVs = sum(prevalence == 1),
    ASVs_prevalence_ge_90pct = sum(prevalence >= 0.90),
    ASVs_prevalence_ge_80pct = sum(prevalence >= 0.80),
    ASVs_prevalence_ge_50pct = sum(prevalence >= 0.50),
    .groups = "drop"
  )

core_summary_by_host

# compares to the higher order zeta
high_order_zeta <- zeta_summary_host %>%
  dplyr::group_by(Predator) %>%
  dplyr::filter(order == max(order)) %>%
  dplyr::select(Predator, order, zeta)

core_summary_by_host %>%
  dplyr::left_join(high_order_zeta, by = "Predator")


## readable labels to match those in differential abundance plots
tax_df_full <- as.data.frame(
  phyloseq::tax_table(ps.major)
) %>%
  tibble::rownames_to_column("category") %>%
  dplyr::mutate(
    category = as.character(category),
    
    ASV_short = dplyr::case_when(
      grepl("^ASV[0-9]+$", category) ~ category,
      TRUE ~ substr(category, 1, 12)
    ),
    
    cat_small = dplyr::case_when(
      !is.na(Genus) &
        Genus != "" &
        !tolower(Genus) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("g__ ", Genus, " | ", ASV_short),
      
      !is.na(Family) &
        Family != "" &
        !tolower(Family) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("f__ ", Family, " | ", ASV_short),
      
      !is.na(Order) &
        Order != "" &
        !tolower(Order) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("o__ ", Order, " | ", ASV_short),
      
      !is.na(Class) &
        Class != "" &
        !tolower(Class) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("c__ ", Class, " | ", ASV_short),
      
      !is.na(Phylum) &
        Phylum != "" &
        !tolower(Phylum) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("p__ ", Phylum, " | ", ASV_short),
      
      TRUE ~ paste0("Unclassified | ", ASV_short)
    )
  )


tax_df_for_core_function <- tax_df_full %>%
  dplyr::select(-category, -ASV_short, -cat_small) %>%
  as.data.frame()

rownames(tax_df_for_core_function) <- tax_df_full$category

# compares
core_from_zeta_groups <- purrr::imap_dfr(
  sample_ids_by_predator,
  function(sample_ids, predator_name) {
    
    make_core_table(
      count_matrix = asv_mat[sample_ids, , drop = FALSE],
      tax_table_df = tax_df_for_core_function,
      prevalence_threshold = 0.80,
      detection_threshold = 0.001
    ) %>%
      dplyr::mutate(
        Predator = predator_name,
        .before = 1
      )
  }
) %>%
  dplyr::left_join(
    tax_df_full %>%
      dplyr::select(category, ASV_short, cat_small),
    by = "category"
  )


core_lookup <- core_from_zeta_groups %>%
  dplyr::select(
    Predator,
    category,
    cat_small,
    n_samples,
    n_present,
    prevalence,
    total_reads,
    mean_rel_abundance_all,
    mean_rel_abundance_present,
    max_rel_abundance,
    core_prevalence,
    core_abundance,
    core_candidate
  ) %>%
  dplyr::mutate(
    core_class = dplyr::case_when(
      prevalence == 1 ~ "Universal core (100%)",
      prevalence >= 0.90 ~ "High-prevalence core (>=90%)",
      prevalence >= 0.80 & core_candidate ~ "Candidate core (>=80%)",
      TRUE ~ "Non-core"
    ),
    core_class = factor(
      core_class,
      levels = c(
        "Universal core (100%)",
        "High-prevalence core (>=90%)",
        "Candidate core (>=80%)",
        "Non-core"
      )
    )
  )


dir.create(
  "Deliverables/ubiome/core_asvs",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  core_lookup,
  "Deliverables/ubiome/core_asvs/ubiome_ASV_prevalence_by_host.csv",
  row.names = FALSE
)

write.csv(
  core_lookup %>%
    dplyr::filter(core_candidate),
  "Deliverables/ubiome/core_asvs/ubiome_candidate_core_ASVs_by_host.csv",
  row.names = FALSE
)

# saves for differential abundance plotting
saveRDS(
  core_lookup,
  file = "./Scripts/ubiome/rdata/ubiome_core_ASV_lookup.rds"
)
