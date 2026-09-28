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
# ??radEmu

# loads in the data
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

# removes negs and mocks
metadata <- metadata %>% drop_na(Predator) 

# optionally removes fin whale, since there is only one sample (two sub-samples of the same sample)
metadata <- metadata %>%
  dplyr::filter(Predator != "fin whale")

asv_mat <- as.data.frame(asv_mat) %>% 
  rownames_to_column(var = "row_id") %>% 
  dplyr::filter(!row_id %in% c("WADE-003-150-A", "WADE-003-150-B",
                               "neg1","mock1-782D", "mock2-783D", "neg2","neg2-2" )) %>% 
  column_to_rownames(var = "row_id")

# # optionally removes spring and winter, since there are few samples
# metadata <- metadata %>%
#   dplyr::filter(season != "Winter")%>%
#   dplyr::filter(season != "Spring")
# asv_mat <- asv_mat %>%
#   rownames_to_column(var = "LabID") %>%
#   dplyr::filter(LabID %in% metadata$LabID) %>%
#   column_to_rownames(var = "LabID")
# 
# # optionally removes south bering from ringed seals (only 2 samples)
# metadata <- metadata %>%
#   dplyr::filter(
#     !(Predator == "ringed seal" & Locale == "South Bering")
#   )
# asv_mat <- asv_mat %>%
#   rownames_to_column(var = "LabID") %>%
#   dplyr::filter(LabID %in% metadata$LabID) %>%
#   column_to_rownames(var = "LabID")

# Samples in metadata but not in asv mat
setdiff(rownames(metadata), rownames(asv_mat))

# Samples in asv mat but not in metadata
setdiff(rownames(asv_mat), rownames(metadata))

# Keep in ps.norep only samples retained in the filtered metadata data frame
ps.pred <- phyloseq::prune_samples(
  phyloseq::sample_names(ps.norep) %in% rownames(metadata),
  ps.norep
)

# Samples retained in ps.norep but absent from ps.pred
samples_only_in_norep <- setdiff(
  phyloseq::sample_names(ps.norep),
  phyloseq::sample_names(ps.pred)
)

samples_only_in_norep

removed_sample_info <- data.frame(
  SampleID = samples_only_in_norep,
  total_counts = phyloseq::sample_sums(ps.norep)[samples_only_in_norep],
  sample_data(ps.norep)[samples_only_in_norep, ],
  check.names = FALSE
)

removed_sample_info %>%
  dplyr::arrange(total_counts)

# # formatting metadata for age class (only bearded and ringed seals)
# metadata_age <- metadata %>%
#   dplyr::filter(Predator %in% c("ringed seal", "bearded seal")) %>%
#   dplyr::mutate(
#     Age_Group = dplyr::case_when(
#       Age_Group %in% c("", "Pending") ~ "Unknown",
#       TRUE ~ Age_Group
#     )
#   )
# 
# metadata_age <- metadata_age %>%
#   dplyr::filter(!Age_Group == "Unknown")
# 
# # formatting asv_mat for age class (only bearded and ringed seals)
# asv_mat_age <- asv_mat %>%
#   rownames_to_column(var = "LabID") %>%
#   dplyr::filter(LabID %in% metadata_age$LabID) %>%
#   column_to_rownames(var = "LabID")

# ------------------------------------------------------------------------------
# FORMATS PS OBJECT TO REMOVE SAMPLES WITH 0 SUMS IN THE OTU TABLE AND NAS IN RELEVANT COVARIATE COLUMNS
# ------------------------------------------------------------------------------
# Inspect sample counts before filtering
phyloseq::nsamples(ps.pred)

# finds samples that have zero abundance
zero_count_pred <- phyloseq::sample_names(ps.pred)[
  phyloseq::sample_sums(ps.pred) == 0
]

zero_count_pred

# Removes samples in zero_count_pred
# Remove samples with zero counts across all ASVs
ps.pred <- phyloseq::prune_samples(
  phyloseq::sample_sums(ps.pred) > 0,
  ps.pred
)

# Now remove ASVs with zero counts across the retained samples
ps.pred <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps.pred) > 0,
  ps.pred
)

# Confirm outcome
phyloseq::nsamples(ps.pred)
phyloseq::ntaxa(ps.pred)

# ------------------------------------------------------------------------------
# CREATE READABLE ASV LABELS FOR radEmu PLOTS
# ------------------------------------------------------------------------------

# Preserve the actual ASV IDs used by the phyloseq object / radEmu model.
# Do NOT do: colnames(asv_mat) <- tax_table(ps.major)$Class
# Classes repeat across multiple ASVs and are not unique identifiers.

tax_df <- as.data.frame(phyloseq::tax_table(ps.pred)) %>%
  tibble::rownames_to_column("category") %>%
  dplyr::mutate(
    category = as.character(category),
    
    ASV_short = dplyr::case_when(
      grepl("^ASV[0-9]+$", category) ~ category,
      TRUE ~ substr(category, 1, 12)
    ),
    
    cat_small = dplyr::case_when(
      !is.na(Species) &
        Species != "" &
        !tolower(Species) %in% c(
          "na", "unknown", "unclassified", "uncultured"
        ) ~ paste0("s__ ", Species, " | ", ASV_short),
      
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


# ---------------------------------------------------------------------------------------
# MODELING
# ---------------------------------------------------------------------------------------

# FITTING A MODEL FOR season---------------------------------------------------------------

# sets the season groups as factors
## beluga whale is the reference group in this case (change the order in levels to  change this)
phyloseq::sample_data(ps.pred)$Predator <- factor(phyloseq::sample_data(ps.pred)$Predator, levels = c("beluga whale","bearded seal", "ringed seal"))

# Fit the seasonal radEmu model once
ch_fit <- emuFit(
  formula = ~ Predator,
  Y = ps.pred,
  run_score_tests = FALSE
)

# Identify which model-matrix columns are host contrasts
host_meta <- data.frame(phyloseq::sample_data(ps.norep))

X_host <- model.matrix(
  ~ Predator,
  data = host_meta
)

host_k <- grep("^Predator", colnames(X_host))

# Specify every season-coefficient × ASV test
test_kj_host <- expand.grid(
  k = host_k,
  j = seq_len(phyloseq::ntaxa(ps.norep))
)

# Run score tests using the model fit
test_all_shost <- emuFit(
  formula = ~ host,
  Y = ps.norep,
  B = ch_fit$B,
  run_score_tests = TRUE,
  test_kj = test_kj_host
)

test_all_host

# plots
p1 <- plot(
  ch_fit,
  taxon_names = tax_df,
  display_taxon_names = TRUE
)$plots

p1

p2 <- plot(
  test_all_season,
  taxon_names = tax_df, 
  display_taxon_names = TRUE)$plots
p2

# SAVES PLOTS--------------------------------------------------------------------------
library(patchwork)

# Generate the radEmu plot list
p2 <- plot(
  test_all_season,
  taxon_names = tax_df,
  display_taxon_names = TRUE
)$plots

# Check the plot names and confirm there are three
names(p2)
length(p2)

# Combine all plots into one horizontal figure
season_plot_panels <- patchwork::wrap_plots(
  p2,
  ncol = 1
)

# View in RStudio
season_plot_panels <- season_plot_panels&
  theme(
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 25),
    plot.subtitle = element_text(size = 21),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 21)
  )

# Save as a high-resolution presentation/publication figure
ggsave(
  filename = "Deliverables/ALL/differential_abundance/pod_radEmu_panels.png",
  plot = season_plot_panels,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)




# FITTING A MODEL FOR POD---------------------------------------------------------------

# sets the season groups as factors
## beluga whale is the reference group in this case (change the order in levels to  change this)
phyloseq::sample_data(ps.16s.major)$pod_from_id<- factor(phyloseq::sample_data(ps.16s.major)$pod_from_id, levels = c("J","K", "L"))

# Fit the seasonal radEmu model once
ch_fit_pod <- emuFit(
  formula = ~ pod_from_id,
  Y = ps.16s.major,
  run_score_tests = FALSE
)

# Identify which model-matrix columns are season contrasts
pod_meta <- data.frame(phyloseq::sample_data(ps.16s.major))

X_pod <- model.matrix(
  ~ pod_from_id,
  data = pod_meta
)

pod_k <- grep("^pod_from_id", colnames(X_pod))

# Specify every season-coefficient × ASV test
test_kj_pod <- expand.grid(
  k = pod_k,
  j = seq_len(phyloseq::ntaxa(ps.16s.major))
)

# Run score tests using the model fit
test_all_pod <- emuFit(
  formula = ~ pod_from_id,
  Y = ps.16s.major,
  B = ch_fit_pod$B,
  run_score_tests = TRUE,
  test_kj = test_kj_pod
)

test_all_pod

# plots
p1 <- plot(
  ch_fit_pod,
  taxon_names = tax_df,
  display_taxon_names = TRUE
)$plots

p1

p2 <- plot(
  test_all_pod,
  taxon_names = tax_df, 
  display_taxon_names = TRUE)$plots
p2

# SAVES PLOTS--------------------------------------------------------------------------
library(patchwork)

# Generate the radEmu plot list
p2 <- plot(
  test_all_pod,
  taxon_names = tax_df,
  display_taxon_names = TRUE
)$plots

# Check the plot names and confirm there are three
names(p2)
length(p2)

# Combine all plots into one horizontal figure
pod_plot_panels <- patchwork::wrap_plots(
  p2,
  ncol = 1
)

# View in RStudio
pod_plot_panels <- pod_plot_panels&
  theme(
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 25),
    plot.subtitle = element_text(size = 21),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 21)
  )

# Save as a high-resolution presentation/publication figure
ggsave(
  filename = "Deliverables/ALL/differential_abundance/pod_radEmu_panels.png",
  plot = pod_plot_panels,
  width = 20,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)

