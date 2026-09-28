# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS DIFFERENTIAL ABUNDANCE METRICS ON UBIOME DATA
# THIS IS THE SECOND STATISTICAL INVESTIGATION 
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
# ------------------------------------------------------------------------------
## radEmu: Using Relative Abundance Data to Estimate of Multiplicative Differences in Mean Absolute Abundance

# sets up the environment
library(radEmu) 
library(phyloseq)
library(tidyverse)
library(patchwork)

# pulls up the tutorial vignettes for radEmu
??radEmu

# loads in the data
load("./Scripts/ubiome/rdata/ubiome.ps-major.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.major), "matrix")

# formats for age group

# reassigns a phyloseq with just ringed and bearded seals
ps.age <- phyloseq::subset_samples(
  ps.major,
  Predator %in% c("ringed seal", "bearded seal")
)

# ensures no 0's remain
ps.age <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps.age) > 0,
  ps.age
)

# recodes 
phyloseq::sample_data(ps.age)$Age_Group <- dplyr::case_when(
  is.na(phyloseq::sample_data(ps.age)$Age_Group) |
    phyloseq::sample_data(ps.age)$Age_Group %in% c("", "Pending") ~ "Unknown",
  TRUE ~ as.character(phyloseq::sample_data(ps.age)$Age_Group)
)

# ------------------------------------------------------------------------------
# CREATE READABLE ASV LABELS FOR radEmu PLOTS
# ------------------------------------------------------------------------------

# Preserve the actual ASV IDs used by the phyloseq object / radEmu model.
# Do NOT do: colnames(asv_mat) <- tax_table(ps.major)$Class
# Classes repeat across multiple ASVs and are not unique identifiers.

tax_df <- as.data.frame(phyloseq::tax_table(ps.major)) %>%
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
        Phylum != "" &
        !tolower(Phylum) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("p__ ", Phylum, " | ", ASV_short),
      
      TRUE ~ paste0("Unclassified | ", ASV_short)
    )
  ) %>%
  dplyr::select(category, cat_small)


# ==============================================================================
# MODELING: radEmu DIFFERENTIAL ABUNDANCE
#
# Analyses:
#   1. Host: all hosts, ~ Predator
#   2. Season within each host: ~ season
#   3. Locale within beluga and ringed seals: ~ Locale
#   4. Age group within each seal host: ~ Age_Group
#
# Important:
#   - Each within-host model is fitted to a separate phyloseq object.
#   - Taxa with zero total abundance are pruned after every sample subset.
#   - Missing/blank grouping values are removed before model fitting.
#   - "Unknown" age values are excluded from the primary biological age analysis.
# ==============================================================================

# ------------------------------------------------------------------------------
# HELPER FUNCTION: SUBSET A PHYLOSEQ OBJECT TO ONE HOST
# ------------------------------------------------------------------------------

make_host_ps <- function(
    ps,
    host_name,
    grouping_var = NULL,
    allowed_levels = NULL
) {
  
  # Extract metadata as a conventional data frame.
  sample_df <- as.data.frame(
    phyloseq::sample_data(ps)
  )
  
  # Keep samples belonging to the requested host.
  keep_host <- rownames(sample_df)[
    as.character(sample_df$Predator) == host_name
  ]
  
  # Use prune_samples(), which accepts evaluated sample IDs.
  ps_sub <- phyloseq::prune_samples(
    keep_host,
    ps
  )
  
  # If requested, remove samples with missing/blank grouping values
  # and optionally retain only specified biological levels.
  if (!is.null(grouping_var)) {
    
    sample_df_sub <- as.data.frame(
      phyloseq::sample_data(ps_sub)
    )
    
    keep_samples <- rownames(sample_df_sub)[
      !is.na(sample_df_sub[[grouping_var]]) &
        trimws(as.character(sample_df_sub[[grouping_var]])) != ""
    ]
    
    if (!is.null(allowed_levels)) {
      keep_samples <- rownames(sample_df_sub)[
        !is.na(sample_df_sub[[grouping_var]]) &
          trimws(as.character(sample_df_sub[[grouping_var]])) != "" &
          as.character(sample_df_sub[[grouping_var]]) %in% allowed_levels
      ]
    }
    
    ps_sub <- phyloseq::prune_samples(
      keep_samples,
      ps_sub
    )
  }
  
  # Remove ASVs absent from every retained sample.
  ps_sub <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(ps_sub) > 0,
    ps_sub
  )
  
  return(ps_sub)
}

# ------------------------------------------------------------------------------
# HELPER FUNCTION: PREPARE A TAXON LABEL TABLE FOR A SUBSETTED PHYLOSEQ OBJECT
# ------------------------------------------------------------------------------

make_tax_df_subset <- function(ps_sub, tax_df_full) {
  
  tax_df_full %>%
    dplyr::filter(
      category %in% phyloseq::taxa_names(ps_sub)
    )
}

# ------------------------------------------------------------------------------
# HELPER FUNCTION: CHECK SAMPLE AND TAXON COUNTS BEFORE MODEL FITTING
# ------------------------------------------------------------------------------

check_ps_for_emu <- function(ps_sub, grouping_var) {
  
  sample_df <- as.data.frame(
    phyloseq::sample_data(ps_sub)
  )
  
  cat("\n--------------------------------------------------\n")
  cat("Grouping variable:", grouping_var, "\n")
  cat("Number of samples:", phyloseq::nsamples(ps_sub), "\n")
  cat("Number of ASVs:", phyloseq::ntaxa(ps_sub), "\n")
  cat("Samples per group:\n")
  print(table(sample_df[[grouping_var]], useNA = "ifany"))
  
  cat("\nSamples with zero reads:\n")
  print(sum(phyloseq::sample_sums(ps_sub) == 0))
  
  cat("\nASVs with zero reads:\n")
  print(sum(phyloseq::taxa_sums(ps_sub) == 0))
  
  cat("--------------------------------------------------\n")
}

# ------------------------------------------------------------------------------
# HELPER FUNCTION: CONFIRM MODEL ASV LABELS MATCH TAXONOMY LOOKUP
# ------------------------------------------------------------------------------

check_taxon_labels <- function(fit, tax_df_subset) {
  
  unmatched_taxa <- fit$coef %>%
    dplyr::distinct(category) %>%
    dplyr::anti_join(
      tax_df_subset,
      by = "category"
    )
  
  if (nrow(unmatched_taxa) == 0) {
    message("Taxonomy-label check passed: all fitted ASVs have readable labels.")
  } else {
    warning(
      "Some fitted ASVs have no label in tax_df_subset."
    )
    print(unmatched_taxa)
  }
}

# ==============================================================================
# 1. DIFFERENTIAL ABUNDANCE ACROSS HOSTS
# ==============================================================================

# All-host object: remove taxa absent across the full data set, if any.
ps.host <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps.major) > 0,
  ps.major
)

# Beluga whale is the reference host.
phyloseq::sample_data(ps.host)$Predator <- factor(
  phyloseq::sample_data(ps.host)$Predator,
  levels = c(
    "beluga whale",
    "bearded seal",
    "ringed seal"
  )
)

tax_df_host <- make_tax_df_subset(
  ps.host,
  tax_df
)

check_ps_for_emu(
  ps.host,
  grouping_var = "Predator"
)

host_fit <- radEmu::emuFit(
  formula = ~ Predator,
  Y = ps.host,
  test_kj = expand.grid(
    k = c(2, 3),
    j = seq_len(phyloseq::ntaxa(ps.host))
  )
)

host_fit

check_taxon_labels(
  host_fit,
  tax_df_host
)

p_host <- plot(
  host_fit,
  taxon_names = tax_df_host,
  display_taxon_names = TRUE
)$plots

p_host

# ==============================================================================
# 2. SEASONAL DIFFERENTIAL ABUNDANCE WITHIN EACH HOST
# ==============================================================================

# Reference level: Winter.
# Note: a comparison is estimable only if the host has samples in Winter and in
# the comparison season. Sparse groups should be interpreted cautiously.

season_levels <- c(
  "Summer",
  "Autumn"
)

# ------------------------------------------------------------------------------
# 2A. BELUGA WHALE: SEASON
# ------------------------------------------------------------------------------

ps.beluga.season <- make_host_ps(
  ps = ps.major,
  host_name = "beluga whale",
  grouping_var = "season",
  allowed_levels = season_levels
)

phyloseq::sample_data(ps.beluga.season)$season <- factor(
  phyloseq::sample_data(ps.beluga.season)$season,
  levels = season_levels
)

tax_df_beluga_season <- make_tax_df_subset(
  ps.beluga.season,
  tax_df
)

check_ps_for_emu(
  ps.beluga.season,
  grouping_var = "season"
)

beluga_season_fit <- radEmu::emuFit(
  formula = ~ season,
  Y = ps.beluga.season,
  test_kj = data.frame(
    k = 2,
    j = seq_len(phyloseq::ntaxa(ps.beluga.season))
  )
)

beluga_season_fit

check_taxon_labels(
  beluga_season_fit,
  tax_df_beluga_season
)

p_beluga_season <- plot(
  beluga_season_fit,
  taxon_names = tax_df_beluga_season,
  display_taxon_names = TRUE,
  title = "Beluga whale: Autumn vs Summer"
)$plots

p_beluga_season

# ------------------------------------------------------------------------------
# 2B. RINGED SEAL: SEASON
# ------------------------------------------------------------------------------

ps.ringed.season <- make_host_ps(
  ps = ps.major,
  host_name = "ringed seal",
  grouping_var = "season",
  allowed_levels = season_levels
)

phyloseq::sample_data(ps.ringed.season)$season <- factor(
  phyloseq::sample_data(ps.ringed.season)$season,
  levels = season_levels
)

tax_df_ringed_season <- make_tax_df_subset(
  ps.ringed.season,
  tax_df
)

check_ps_for_emu(
  ps.ringed.season,
  grouping_var = "season"
)

ringed_season_fit <- radEmu::emuFit(
  formula = ~ season,
  Y = ps.ringed.season,
  test_kj = data.frame(
    k = 2,
    j = seq_len(phyloseq::ntaxa(ps.ringed.season))
  )
)

ringed_season_fit

check_taxon_labels(
  ringed_season_fit,
  tax_df_ringed_season
)

p_ringed_season <- plot(
  ringed_season_fit,
  taxon_names = tax_df_ringed_season,
  display_taxon_names = TRUE,
  title = "Ringed seal: Autumn vs Summer"
)$plots

p_ringed_season

# ------------------------------------------------------------------------------
# 2C. BEARDED SEAL: SEASON
# ------------------------------------------------------------------------------

ps.bearded.season <- make_host_ps(
  ps = ps.major,
  host_name = "bearded seal",
  grouping_var = "season",
  allowed_levels = season_levels
)

phyloseq::sample_data(ps.bearded.season)$season <- factor(
  phyloseq::sample_data(ps.bearded.season)$season,
  levels = season_levels
)

tax_df_bearded_season <- make_tax_df_subset(
  ps.bearded.season,
  tax_df
)

check_ps_for_emu(
  ps.bearded.season,
  grouping_var = "season"
)

bearded_season_fit <- radEmu::emuFit(
  formula = ~ season,
  Y = ps.bearded.season,
  test_kj = data.frame(
    k = 2,
    j = seq_len(phyloseq::ntaxa(ps.bearded.season))
  )
)

bearded_season_fit

check_taxon_labels(
  bearded_season_fit,
  tax_df_bearded_season
)

p_bearded_season <- plot(
  bearded_season_fit,
  taxon_names = tax_df_bearded_season,
  display_taxon_names = TRUE,
  title = "Bearded seal: Autumn vs Summer"
)$plots

p_bearded_season

# ------------------------------------------------------------------------------
# SEASON PLOTS: COMBINE BY HOST
# ------------------------------------------------------------------------------

# fixed axes labels 
p_bearded_season[[1]] <- p_bearded_season[[1]] +
  labs(x = "Estimated fold difference in mean ASV abundance", y = "ASV")

p_beluga_season[[1]] <- p_beluga_season[[1]] +
  labs(x = "Estimated fold difference in mean ASV abundance", y = "ASV")

p_ringed_season[[1]] <- p_ringed_season[[1]] +
  labs(x = "Estimated fold difference in mean ASV abundance", y = "ASV")

season_plots_by_host <- (
  p_bearded_season[[1]] |
    p_beluga_season[[1]] |
    p_ringed_season[[1]]
) +
  patchwork::plot_annotation(
    title = "Seasonal differences in relative abundance by host"
  ) &
  theme(
    text = element_text(size = 12),
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(t = 8, b = 12)
    ), legend.position = "none"
  )

season_plots_by_host
   
ggsave(
  filename = "Deliverables/ubiome/diffabund/DIFFABUND_ubiome_by-season_each-host.png",
  plot = season_plots_by_host,
  width = 30,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ==============================================================================
# 3. LOCALE DIFFERENTIAL ABUNDANCE WITHIN HOST
# ==============================================================================

# Analyze only beluga whales and ringed seals.
# Bearded seals are excluded because all samples belong to Arctic locale, so
# there is no within-host locale contrast.

locale_levels <- c(
  "Arctic",
  "South Bering",
  "Cook Inlet"
)

# ------------------------------------------------------------------------------
# 3A. BELUGA WHALE: LOCALE
# ------------------------------------------------------------------------------

ps.beluga.locale <- make_host_ps(
  ps = ps.major,
  host_name = "beluga whale",
  grouping_var = "Locale",
  allowed_levels = locale_levels
)

phyloseq::sample_data(ps.beluga.locale)$Locale <- factor(
  phyloseq::sample_data(ps.beluga.locale)$Locale,
  levels = locale_levels
)

# Drop unused factor levels. This is particularly important because the
# host-specific data may not contain every global locale category.
phyloseq::sample_data(ps.beluga.locale)$Locale <- droplevels(
  phyloseq::sample_data(ps.beluga.locale)$Locale
)

tax_df_beluga_locale <- make_tax_df_subset(
  ps.beluga.locale,
  tax_df
)

check_ps_for_emu(
  ps.beluga.locale,
  grouping_var = "Locale"
)

beluga_locale_fit <- radEmu::emuFit(
  formula = ~ Locale,
  Y = ps.beluga.locale,
  test_kj = expand.grid(
    k = c(2, 3),
    j = seq_len(phyloseq::ntaxa(ps.beluga.locale))
  )
)

beluga_locale_fit

check_taxon_labels(
  beluga_locale_fit,
  tax_df_beluga_locale
)

p_beluga_locale <- plot(
  beluga_locale_fit,
  taxon_names = tax_df_beluga_locale,
  display_taxon_names = TRUE
)$plots

p_beluga_locale

p_beluga_locale[[1]] <- p_beluga_locale[[1]] +
  labs(
    title = "Beluga whale: South Bering vs Arctic",
    x = "Estimated fold difference in mean ASV abundance",
    y = "ASV"
  )

p_beluga_locale[[2]] <- p_beluga_locale[[2]] +
  labs(
    title = "Beluga whale: Cook Inlet vs Arctic",
    x = "Estimated fold difference in mean ASV abundance",
    y = "ASV"
  )

p_beluga_locale

# ------------------------------------------------------------------------------
# 3B. RINGED SEAL: LOCALE
## can't do a comparison for ringed or bearded seals 
## because of low sample in South Bering and Cook Inlet for both
# ------------------------------------------------------------------------------

# Ringed seal data may only contain Arctic and South Bering. After droplevels(),
# Arctic remains the reference and South Bering is the single comparison.

# ps.ringed.locale <- make_host_ps(
#   ps = ps.major,
#   host_name = "ringed seal",
#   grouping_var = "Locale",
#   allowed_levels = locale_levels
# )
# 
# phyloseq::sample_data(ps.ringed.locale)$Locale <- factor(
#   phyloseq::sample_data(ps.ringed.locale)$Locale,
#   levels = locale_levels
# )
# 
# phyloseq::sample_data(ps.ringed.locale)$Locale <- droplevels(
#   phyloseq::sample_data(ps.ringed.locale)$Locale
# )
# 
# tax_df_ringed_locale <- make_tax_df_subset(
#   ps.ringed.locale,
#   tax_df
# )
# 
# check_ps_for_emu(
#   ps.ringed.locale,
#   grouping_var = "Locale"
# )
# 
# ringed_locale_fit <- radEmu::emuFit(
#   formula = ~ Locale,
#   Y = ps.ringed.locale,
#   test_kj = data.frame(
#     k = 2,
#     j = 1
#   )
# )
# 
# ringed_locale_fit
# 
# check_taxon_labels(
#   ringed_locale_fit,
#   tax_df_ringed_locale
# )
# 
# p_ringed_locale <- plot(
#   ringed_locale_fit,
#   taxon_names = tax_df_ringed_locale,
#   display_taxon_names = TRUE
#   
# )$plots
# 
# p_ringed_locale
# 
# p_ringed_locale[[1]] <- p_ringed_locale[[1]] +
#   labs(
#     title = "Ringed seal: South Bering vs Arctic",
#     x = "Estimated fold difference in mean ASV abundance",
#     y = "ASV"
#   )

# ------------------------------------------------------------------------------
# LOCALE PLOTS: COMBINE BY HOST
# ------------------------------------------------------------------------------

beluga_locale_plots <- patchwork::wrap_plots(
  p_beluga_locale,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Differential abundance by locale within beluga whales",
    tag_levels = "A"
  )

# ringed_locale_plots <- patchwork::wrap_plots(
#   p_ringed_locale,
#   ncol = 1
# ) +
#   patchwork::plot_annotation(
#     title = "Differential abundance by locale within ringed seals",
#     tag_levels = "A"
#   )

locale_plots_by_host <- (
 beluga_locale_plots #/
  #   ringed_locale_plots
)

locale_plots_by_host

ggsave(
  "Deliverables/ubiome/diffabund/DIFFABUND_ubiome_by-locale_each-host.png",
  plot = locale_plots_by_host,
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

# ==============================================================================
# 4. AGE-GROUP DIFFERENTIAL ABUNDANCE WITHIN EACH SEAL HOST
# ==============================================================================

# Primary age comparison:
#   Pup vs Non-pup
#
# Unknown / pending age data are omitted so that the model answers the
# biological age-group question rather than a known-vs-missing-age comparison.

age_levels <- c(
  "Non-pup",
  "Pup"
)

# ------------------------------------------------------------------------------
# 4A. RINGED SEAL: AGE GROUP
# ------------------------------------------------------------------------------

ps.ringed.age <- make_host_ps(
  ps = ps.major,
  host_name = "ringed seal",
  grouping_var = "Age_Group",
  allowed_levels = age_levels
)

phyloseq::sample_data(ps.ringed.age)$Age_Group <- factor(
  phyloseq::sample_data(ps.ringed.age)$Age_Group,
  levels = age_levels
)

tax_df_ringed_age <- make_tax_df_subset(
  ps.ringed.age,
  tax_df
)

check_ps_for_emu(
  ps.ringed.age,
  grouping_var = "Age_Group"
)

ringed_age_fit <- radEmu::emuFit(
  formula = ~ Age_Group,
  Y = ps.ringed.age,
  test_kj = data.frame(
    k = 2,
    j = seq_len(phyloseq::ntaxa(ps.ringed.age))
  )
)

ringed_age_fit

check_taxon_labels(
  ringed_age_fit,
  tax_df_ringed_age
)

p_ringed_age <- plot(
  ringed_age_fit,
  taxon_names = tax_df_ringed_age,
  display_taxon_names = TRUE
)$plots

p_ringed_age

p_ringed_age[[1]] <- p_ringed_age[[1]] +
  labs(
    title = "Ringed seal: Pup vs Non-Pup",
    x = "Estimated fold difference in mean ASV abundance",
    y = "ASV"
  )

# ------------------------------------------------------------------------------
# 4B. BEARDED SEAL: AGE GROUP
# ------------------------------------------------------------------------------

ps.bearded.age <- make_host_ps(
  ps = ps.major,
  host_name = "bearded seal",
  grouping_var = "Age_Group",
  allowed_levels = age_levels
)

phyloseq::sample_data(ps.bearded.age)$Age_Group <- factor(
  phyloseq::sample_data(ps.bearded.age)$Age_Group,
  levels = age_levels
)

tax_df_bearded_age <- make_tax_df_subset(
  ps.bearded.age,
  tax_df
)

check_ps_for_emu(
  ps.bearded.age,
  grouping_var = "Age_Group"
)

bearded_age_fit <- radEmu::emuFit(
  formula = ~ Age_Group,
  Y = ps.bearded.age,
  test_kj = data.frame(
    k = 2,
    j = seq_len(phyloseq::ntaxa(ps.bearded.age))
  )
)

bearded_age_fit

check_taxon_labels(
  bearded_age_fit,
  tax_df_bearded_age
)

p_bearded_age <- plot(
  bearded_age_fit,
  taxon_names = tax_df_bearded_age,
  display_taxon_names = TRUE
)$plots

p_bearded_age

p_bearded_age[[1]] <- p_bearded_age[[1]] +
  labs(
    title = "Bearded seal: Pup vs Non-Pup",
    x = "Estimated fold difference in mean ASV abundance",
    y = "ASV"
  )


# ------------------------------------------------------------------------------
# AGE-GROUP PLOTS: COMBINE BY HOST
# ------------------------------------------------------------------------------

ringed_age_plots <- patchwork::wrap_plots(
  p_ringed_age,
  ncol = 1
) +
  patchwork::plot_annotation(
    title = "Differential abundance by age group within ringed seals",
    tag_levels = "A"
  )

bearded_age_plots <- patchwork::wrap_plots(
  p_bearded_age,
  ncol = 1
) +
  patchwork::plot_annotation(
    title = "Differential abundance by age group within bearded seals",
    tag_levels = "A"
  )

age_plots_by_host <- (
  ringed_age_plots |
    bearded_age_plots
)

age_plots_by_host

ggsave(
  "Deliverables/ubiome/diffabund/DIFFABUND_ubiome_by-age_each-seal-host.png",
  plot = age_plots_by_host,
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

# ==============================================================================
# OPTIONAL: SAVE ALL WITHIN-HOST RESULTS IN ONE LARGE FIGURE
# ==============================================================================

all_within_host_diffabund_plots <- (
  season_plots_by_host /
    locale_plots_by_host /
    age_plots_by_host
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

all_within_host_diffabund_plots

ggsave(
  "Deliverables/ubiome/diffabund/DIFFABUND_ubiome_within-host_all.png",
  plot = all_within_host_diffabund_plots,
  width = 24,
  height = 38,
  units = "in",
  dpi = 300
)

# ==============================================================================
# END WITHIN-HOST radEmu DIFFERENTIAL-ABUNDANCE SECTION
# ==============================================================================

# ==============================================================================
# BEGIN, CORE MICROBES SECTION
# ZETADIV MUST BE RUN BEFORE THIS; CORE MB IS BASED ON ZETADIV
# ==============================================================================


core_lookup <- readRDS("./Scripts/ubiome/rdata/ubiome_core_ASV_lookup.rds")

# adds in columns to note which host an ASV is "core" in
core_any_host <- core_lookup %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(
    cat_small = dplyr::first(cat_small),
    core_in_any_host = any(core_candidate),
    universal_in_any_host = any(core_class == "Universal core (100%)"),
    highest_core_class = dplyr::case_when(
      any(core_class == "Universal core (100%)") ~ "Universal core in ≥1 host",
      any(core_class == "High-prevalence core (>=90%)") ~ "High-prevalence core in ≥1 host",
      any(core_class == "Candidate core (>=80%)") ~ "Candidate core in ≥1 host",
      TRUE ~ "Non-core"
    ),
    hosts_core = dplyr::if_else(
      any(core_candidate),
      paste(
        unique(as.character(Predator[core_candidate])),
        collapse = "; "
      ),
      NA_character_
    ),
    max_prevalence = max(prevalence),
    .groups = "drop"
  )

# adds core labels to coefficients
add_core_status <- function(
    fit,
    tax_labels,
    core_lookup,
    host_name = NULL
) {
  
  coef_df <- as.data.frame(fit$coef) %>%
    dplyr::mutate(category = as.character(category)) %>%
    dplyr::left_join(
      tax_labels %>%
        dplyr::mutate(category = as.character(category)) %>%
        dplyr::select(category, cat_small),
      by = "category"
    )
  
  if (is.null(host_name)) {
    
    # Make one row per ASV, recording every host where it is a candidate core.
    core_host_lookup <- core_lookup %>%
      dplyr::filter(core_candidate) %>%
      dplyr::mutate(
        Predator = as.character(Predator)
      ) %>%
      dplyr::group_by(category) %>%
      dplyr::summarise(
        core_hosts = paste(
          sort(unique(Predator)),
          collapse = "; "
        ),
        n_core_hosts = dplyr::n_distinct(Predator),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        core_host_shape = dplyr::case_when(
          n_core_hosts == 3 ~ "All hosts",
          n_core_hosts > 1 ~ "Multiple hosts",
          core_hosts == "bearded seal" ~ "Bearded seal",
          core_hosts == "ringed seal" ~ "Ringed seal",
          core_hosts == "beluga whale" ~ "Beluga whale",
          TRUE ~ "Non-core"
        )
      )
    
    coef_df <- coef_df %>%
      dplyr::left_join(
        core_any_host %>%
          dplyr::select(
            category,
            core_in_any_host,
            universal_in_any_host,
            highest_core_class,
            hosts_core,
            max_prevalence
          ),
        by = "category"
      ) %>%
      dplyr::left_join(
        core_host_lookup %>%
          dplyr::select(
            category,
            core_hosts,
            n_core_hosts,
            core_host_shape
          ),
        by = "category"
      ) %>%
      dplyr::mutate(
        core_status = highest_core_class,
        core_flag = core_in_any_host,
        core_host_shape = dplyr::coalesce(
          core_host_shape,
          "Non-core"
        )
      )
    
  } else {
    
    coef_df <- coef_df %>%
      dplyr::left_join(
        core_lookup %>%
          dplyr::filter(Predator == host_name) %>%
          dplyr::select(
            category,
            prevalence,
            mean_rel_abundance_all,
            core_candidate,
            core_class
          ),
        by = "category"
      ) %>%
      dplyr::mutate(
        core_status = as.character(core_class),
        core_flag = core_candidate,
        
        # For within-host plots, all candidate-core taxa belong to this host.
        core_host_shape = dplyr::if_else(
          core_candidate,
          as.character(host_name),
          "Non-core"
        )
      )
  }
  
  coef_df %>%
    dplyr::mutate(
      cat_small = dplyr::coalesce(cat_small, category),
      core_status = dplyr::coalesce(core_status, "Non-core"),
      core_flag = dplyr::coalesce(core_flag, FALSE),
      core_host_shape = dplyr::coalesce(
        core_host_shape,
        "Non-core"
      ),
      
      # Control both legend order and shape assignment.
      core_host_shape = factor(
        core_host_shape,
        levels = c(
          "Non-core",
          "Bearded seal",
          "Ringed seal",
          "Beluga whale",
          "Multiple hosts",
          "All hosts"
        )
      )
      )
}

# plot 
plot_emu_core <- function(
    coef_df,
    plot_title = NULL,
    max_taxa = 25
) {
  
  plot_df <- coef_df %>%
    dplyr::filter(
      !is.na(estimate),
      !is.na(lower),
      !is.na(upper)
    ) %>%
    dplyr::slice_max(
      order_by = abs(estimate),
      n = max_taxa,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(
      cat_small = forcats::fct_reorder(
        cat_small,
        estimate
      ),
      
      core_status = factor(
        core_status,
        levels = c(
          "Universal core in ≥1 host",
          "High-prevalence core in ≥1 host",
          "Candidate core in ≥1 host",
          "Non-core"
        )
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = estimate,
      y = cat_small,
      colour = core_status,
      shape = core_host_shape
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.8
    ) +
    geom_errorbarh(
      aes(
        xmin = lower,
        xmax = upper
      ),
      height = 0.18,
      linewidth = 1
    ) +
    geom_point(
      size = 4.8,
      stroke = 0.9
    ) +
    scale_colour_manual(
      values = c(
        "Universal core in ≥1 host" = "#542788",
        "High-prevalence core in ≥1 host" = "#2C7FB8",
        "Candidate core in ≥1 host" = "#41AB5D",
        "Non-core" = "grey65"
      ),
      breaks = c(
        "Universal core in ≥1 host",
        "High-prevalence core in ≥1 host",
        "Candidate core in ≥1 host",
        "Non-core"
      ),
      labels = c(
        "Universal core in ≥1 host" = "100% prevalence",
        "High-prevalence core in ≥1 host" = "≥90% prevalence",
        "Candidate core in ≥1 host" = "≥80% prevalence",
        "Non-core" = "Not core"
      ),
      name = "Core prevalence threshold",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "Non-core" = 16,
        "Bearded seal" = 17,
        "Ringed seal" = 15,
        "Beluga whale" = 18,
        "Multiple hosts" = 8,
        "All hosts" = 4
      ),
      breaks = c(
        "Bearded seal",
        "Ringed seal",
        "Beluga whale",
        "Multiple hosts",
        "All hosts"
      ),
      name = "Candidate core host(s)",
      drop = FALSE
    ) +
    labs(
      title = plot_title,
      x = "Estimated log fold-difference in mean absolute abundance",
      y = NULL
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 26,
        hjust = 0.5,
        margin = margin(b = 12)
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 19,
        margin = margin(t = 12)
      ),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold",
        size = 17
      ),
      legend.text = element_text(size = 15),
      legend.box = "vertical",
      legend.spacing.y = unit(0.35, "cm")
    )
}

# applies to host model: 
host_coef_core <- add_core_status(
  fit = host_fit,
  tax_labels = tax_df_host,
  core_lookup = core_lookup,
  host_name = NULL
)

head(host_coef_core)

# plots 
host_core_plots <- host_coef_core %>%
  dplyr::filter(covariate != "(Intercept)") %>%
  dplyr::group_split(covariate) %>%
  purrr::map(
    ~ plot_emu_core(
      coef_df = .x,
      plot_title = unique(.x$covariate),
      max_taxa = 25
    )
  )

host_core_plots

host_core_plot_combined <- patchwork::wrap_plots(
  host_core_plots,
  ncol = 2
) +
  patchwork::plot_annotation(
    title = "Differential abundance across host species",
    subtitle = "Triangles and color indicate host-core ASV status"
  )

host_core_plot_combined

ggsave(
  "Deliverables/ubiome/diffabund/CORE_ubiome_within-host_all.png",
  plot = host_core_plot_combined,
  width = 38,
  height = 24,
  units = "in",
  dpi = 300
)

host_coef_core %>%
  dplyr::filter(core_host_shape == "Multiple hosts") %>%
  dplyr::distinct(cat_small, core_hosts)

# saves environment
getwd()
save.image(file = "./Scripts/ubiome/rdata/radEmu_BYHOST.RData")


host_core_plot_combined <- patchwork::wrap_plots(
  host_core_plots,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Differential abundance across host species",
    subtitle = "Triangles and color indicate host-core ASV status",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 38,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = 24,
        hjust = 0.5,
        margin = margin(b = 22)
      )
    )
  ) &
  theme(
    text = element_text(size = 22),
    
    # Individual panel titles
    plot.title = element_text(
      face = "bold",
      size = 26,
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    
    # Axis text
    axis.title.x = element_text(
      face = "bold",
      size = 20,
      margin = margin(t = 12)
    ),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 16),
    
    # Shared legend
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold",
      size = 19
    ),
    legend.text = element_text(size = 17),
    legend.key.height = unit(0.7, "cm"),
    legend.key.width = unit(0.8, "cm")
  )

host_core_plot_combined
