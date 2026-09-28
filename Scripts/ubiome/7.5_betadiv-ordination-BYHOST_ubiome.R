# ------------------------------------------------------------------------------
# COMPUTES AND INTERPRETS BETA DIVERSITY METRICS ON UBIOME DATA
# THIS IS THE SECOND STATISTICAL INVESTIGATION 
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
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

head(rownames(asv_mat))
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
# Aitchison Distances
# ------------------------------------------------------------------------------
# # STATS : AITCHISON PERMANOVA ON RELATIVE ABUNDANCES 
# # NOT RAREFIED; USED CLR METHOD PER GLOOR ET AL., 2017
# FROM WIKIPEDIA: The null hypothesis states that the centroids (averages) and/or the dispersion (spread) of the groups are equivalent in the multivariate space.
## H0 - there is no difference in the pre samples vs post samples
## HA - there is a difference in the pre samples vs post samples

# Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"
robust_aitch_dist <- vegdist(asv_mat, method = "robust.aitchison")
head(rownames(robust_aitch_dist))

robust_aitch_dist_age <- vegdist(asv_mat_age, method = "robust.aitchison")

# Run the PERMANOVA-------------------------------------------------------------

# BY HOST
permanova_aitch_host <- adonis2(robust_aitch_dist ~ Predator, 
                                data = metadata, 
                                permutations = 999)

print(permanova_aitch_host) # p-value = 0.001

metadata <- data.frame(sample_data(ps.norep)) %>%
  dplyr::filter(LabID %in% labels(robust_aitch_dist)) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist)))

stopifnot(
  identical(metadata$LabID, labels(robust_aitch_dist))
)


pairwise_permanova <- function(dist_obj, metadata, group_var,
                               sample_id = "LabID",
                               permutations = 999,
                               p_adjust = "BH") {
  
  # Ensure metadata and distance labels are identically ordered
  metadata <- metadata %>%
    dplyr::filter(.data[[sample_id]] %in% labels(dist_obj)) %>%
    dplyr::arrange(match(.data[[sample_id]], labels(dist_obj)))
  
  stopifnot(identical(metadata[[sample_id]], labels(dist_obj)))
  
  groups <- sort(unique(metadata[[group_var]]))
  pairs <- combn(groups, 2, simplify = FALSE)
  
  purrr::map_dfr(pairs, function(pair) {
    
    meta_sub <- metadata %>%
      dplyr::filter(.data[[group_var]] %in% pair)
    
    dist_sub <- as.dist(
      as.matrix(dist_obj)[meta_sub[[sample_id]], meta_sub[[sample_id]]]
    )
    
    model <- vegan::adonis2(
      dist_sub ~ group,
      data = data.frame(group = meta_sub[[group_var]]),
      permutations = permutations
    )
    
    tibble::tibble(
      group1 = pair[1],
      group2 = pair[2],
      n1 = sum(meta_sub[[group_var]] == pair[1]),
      n2 = sum(meta_sub[[group_var]] == pair[2]),
      pseudo_F = model$F[1],
      R2 = model$R2[1],
      p = model$`Pr(>F)`[1]
    )
  }) %>%
    dplyr::mutate(
      p_BH = p.adjust(p, method = p_adjust),
      significance = dplyr::case_when(
        p_BH < 0.001 ~ "***",
        p_BH < 0.01 ~ "**",
        p_BH < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::arrange(p_BH, p)
}


set.seed(123)

pairwise_aitch_host <- pairwise_permanova(
  dist_obj = robust_aitch_dist,
  metadata = metadata,
  group_var = "Predator",
  sample_id = "LabID",
  permutations = 9999,
  p_adjust = "BH"
)

pairwise_aitch_host

# tests for dispersion
disp <- vegan::betadisper(
  robust_aitch_dist,
  group = metadata$Predator
)

set.seed(123)
vegan::permutest(disp, permutations = 9999) # p-value = 0.1382
## Multivariate dispersion did not differ significantly among host groups
## (PERMDISP: F = 2.-118, p = 0.0517), although the result was borderline.
## Therefore, the significant PERMANOVA result is not accompanied by formal
## evidence of heterogeneous dispersion at alpha = 0.05.

## However, because this result is so borderline, caution:
## this means that significant results in differences may not be due to 
## multivariate microbial compositional differences,
## but may be due to unequal dispersion

# HOST BOX PLOTS-----------------------------------------------------------------------

# PERMDIST
## this box plot will show within-host beta diversity, 
## or dispersion: whether samples within each predator group are similarly variable.

# creates a usable column from robust aitchison distances
## disp$distances is the distance from each sample to its host-group centroid

disp_plot_df <- metadata %>%
  dplyr::mutate(
    Distance_to_centroid = disp$distances
  )

stopifnot(
  nrow(disp_plot_df) == length(disp$distances),
  identical(as.character(disp_plot_df$Predator), as.character(disp$group))
)
# Larger values mean samples are more dispersed/heterogeneous within that host group
box_plot_host <- ggplot(
  disp_plot_df,
  aes(
    x = Predator,
    y = Distance_to_centroid,
    fill = Predator
  )
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
  labs(
    x = "Host",
    y = "Distance to group centroid",
    title = "Multivariate dispersion by host"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

box_plot_host

# WITHIN VS BETWEEN-HOST DISTANCES
# creates a longform df of robust aitchison distances
pairwise_dist_df <- as.matrix(robust_aitch_dist) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("LabID_1") %>%
  tidyr::pivot_longer(
    cols = -LabID_1,
    names_to = "LabID_2",
    values_to = "Robust_Aitchison"
  ) %>%
  dplyr::filter(LabID_1 < LabID_2) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(LabID, Predator) %>%
      dplyr::rename(
        LabID_1 = LabID,
        Predator_1 = Predator
      ),
    by = "LabID_1"
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(LabID, Predator) %>%
      dplyr::rename(
        LabID_2 = LabID,
        Predator_2 = Predator
      ),
    by = "LabID_2"
  ) %>%
  dplyr::mutate(
    comparison = dplyr::case_when(
      Predator_1 == Predator_2 ~ paste("Within", Predator_1),
      TRUE ~ paste(
        pmin(Predator_1, Predator_2),
        "vs",
        pmax(Predator_1, Predator_2)
      )
    )
  )

# plots
pairwise_dist_df %>%
  ggplot(
    aes(
      x = comparison,
      y = Robust_Aitchison,
      fill = comparison
    )
  ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(
    width = 0.15,
    alpha = 0.15,
    size = 1
  ) +
  labs(
    x = NULL,
    y = "Pairwise robust Aitchison distance",
    title = "Within- and between-host microbiome dissimilarity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "none"
  )

# some are significantly different. 
# bearded seals and beluga whales (BH adjusted p-value = 0.00015)
# ringed seals and beluga whales (BH adjusted p-value = 0.00015)
# bearded seals and ringed seals are not significantly different (BH adjusted p-value = 0.0519)


# Run the PERMANOVA BY SEASON----------------------------------------------------------------------

## cleans metadata table: 
# removes NAs in season
metadata_season <- metadata %>% drop_na(season) 

# cleans avs table based on metadata_seasons
asv_mat_season <- as.data.frame(asv_mat) %>% 
  rownames_to_column(var = "row_id") %>% 
  dplyr::filter(!row_id %in% c("WADE-003-013", "WADE-003-148", "WADE-003-149")) %>% 
  column_to_rownames(var = "row_id")

# separates by host
bearded_metadata_season <- metadata_season %>%
  dplyr::filter(Predator == "bearded seal")
asv_mat_bearded_season <- asv_mat_season %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% bearded_metadata_season$LabID) %>%
  column_to_rownames(var = "LabID")


beluga_metadata_season <- metadata_season %>%
  dplyr::filter(Predator == "beluga whale")
asv_mat_beluga_season <- asv_mat_season %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% beluga_metadata_season$LabID) %>%
  column_to_rownames(var = "LabID")

ringed_metadata_season <- metadata_season %>%
  dplyr::filter(Predator == "ringed seal")
asv_mat_ringed_season <- asv_mat_season %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% ringed_metadata_season$LabID) %>%
  column_to_rownames(var = "LabID")

# RERUN Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"
robust_aitch_dist_bearded_season <- vegdist(asv_mat_bearded_season, method = "robust.aitchison")
head(rownames(robust_aitch_dist_bearded_season))

bearded_permanova_aitch_season <- adonis2(robust_aitch_dist_bearded_season ~ season, 
                                  data = bearded_metadata_season, 
                                  permutations = 999)

print(bearded_permanova_aitch_season) # p-value = 0.074 not significant


robust_aitch_dist_beluga_season <- vegdist(asv_mat_beluga_season, method = "robust.aitchison")
head(rownames(robust_aitch_dist_beluga_season))

beluga_permanova_aitch_season <- adonis2(robust_aitch_dist_beluga_season ~ season, 
                                          data = beluga_metadata_season, 
                                          permutations = 999)

print(beluga_permanova_aitch_season) # p-value = 0.234 not significant

robust_aitch_dist_ringed_season <- vegdist(asv_mat_ringed_season, method = "robust.aitchison")
head(rownames(robust_aitch_dist_ringed_season))

ringed_permanova_aitch_season <- adonis2(robust_aitch_dist_ringed_season ~ season, 
                                         data = ringed_metadata_season, 
                                         permutations = 999)

print(ringed_permanova_aitch_season) # p-value = 0.12 not significant

# 
# 
# metadata_season <- data.frame(sample_data(ps.major)) %>%
#   dplyr::filter(LabID %in% labels(robust_aitch_dist_season)) %>%
#   dplyr::arrange(match(LabID, labels(robust_aitch_dist_season)))
# 
# stopifnot(
#   identical(metadata_season$LabID, labels(robust_aitch_dist_season))
# )
# 
# set.seed(123)
# 
# pairwise_aitch_season <- pairwise_permanova(
#   dist_obj = robust_aitch_dist_season,
#   metadata = metadata_season,
#   group_var = "season",
#   sample_id = "LabID",
#   permutations = 9999,
#   p_adjust = "BH"
# )
# 
# pairwise_aitch_season
# 
# # tests for dispersion
# disp <- vegan::betadisper(
#   robust_aitch_dist_season,
#   group = metadata_season$season
# )
# 
# set.seed(123)
# vegan::permutest(disp, permutations = 9999) # p-value = 0.2476
# ## Multivariate dispersion did not differ significantly among host groups
# 

# No significant differences between seasons


# BY LOCALE----------------------------------------------------------------------------
# removes NAs in season
metadata_Locale <- metadata %>% drop_na(Locale) 

# cleans avs table based on metadata_seasons
asv_mat_Locale <- as.data.frame(asv_mat) %>% 
  rownames_to_column(var = "row_id") %>% 
  # dplyr::filter(!row_id %in% c("WADE-003-013", "WADE-003-148", "WADE-003-149")) %>% 
  column_to_rownames(var = "row_id")


beluga_metadata_Locale <- metadata_Locale %>%
  dplyr::filter(Predator == "beluga whale")
asv_mat_beluga_Locale <- asv_mat_Locale %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% beluga_metadata_Locale$LabID) %>%
  column_to_rownames(var = "LabID")

ringed_metadata_Locale <- metadata_Locale %>%
  dplyr::filter(Predator == "ringed seal")
asv_mat_ringed_Locale <- asv_mat_Locale %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% ringed_metadata_Locale$LabID) %>%
  column_to_rownames(var = "LabID")

# RERUN Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"

#robust_aitch_dist_bearded_Locale <- vegdist(asv_mat_bearded_Locale, method = "robust.aitchison")
## BEARDED SEAL SAMPLES ARE ONLY IN ARCTIC

robust_aitch_dist_beluga_Locale <- vegdist(asv_mat_beluga_Locale, method = "robust.aitchison")
head(rownames(robust_aitch_dist_beluga_Locale))

beluga_permanova_aitch_Locale <- adonis2(robust_aitch_dist_beluga_Locale ~ Locale, 
                                         data = beluga_metadata_Locale, 
                                         permutations = 999)

print(beluga_permanova_aitch_Locale) # p-value = 0.006 significant


# robust_aitch_dist_ringed_Locale <- vegdist(asv_mat_ringed_Locale, method = "robust.aitchison")
# head(rownames(robust_aitch_dist_ringed_Locale))
# 
# ringed_permanova_aitch_Locale <- adonis2(robust_aitch_dist_ringed_Locale ~ Locale, 
#                                          data = ringed_metadata_Locale, 
#                                          permutations = 999)
# 
# print(ringed_permanova_aitch_Locale) # p-value = 0.911 significant
# 


# metadata_beluga_Locale <- data.frame(sample_data(ps.major)) %>%
#   dplyr::filter(LabID %in% labels(robust_aitch_dist_beluga_Locale)) %>%
#   dplyr::arrange(match(LabID, labels(robust_aitch_dist_beluga_Locale)))
# 
# stopifnot(
#   identical(metadata$LabID, labels(robust_aitch_dist))
# )


set.seed(123)

beluga_Locale_pairwise_aitch_locale <- pairwise_permanova(
  dist_obj = robust_aitch_dist_beluga_Locale,
  metadata = beluga_metadata_Locale,
  group_var = "Locale",
  sample_id = "LabID",
  permutations = 9999,
  p_adjust = "BH"
)

beluga_Locale_pairwise_aitch_locale

# tests for dispersion
disp <- vegan::betadisper(
  robust_aitch_dist_beluga_Locale,
  group = beluga_metadata_Locale$Locale
)

set.seed(123)
vegan::permutest(disp, permutations = 9999) # p-value = 0.0278
## Multivariate dispersion did differ significantly among locale groups
## Therefore, the significant PERMANOVA result is accompanied by formal
## evidence of heterogeneous dispersion at alpha = 0.05.


# For belugas, Cook Inlet was significantly different from South Bering and Arctic DB

# BY AGE------------------------------------------------------------------------------

bearded_metadata_age <- metadata_age %>%
  dplyr::filter(Predator == "bearded seal")
asv_mat_bearded_age <- asv_mat_age %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% bearded_metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

ringed_metadata_age <- metadata_age %>%
  dplyr::filter(Predator == "ringed seal")
asv_mat_ringed_age <- asv_mat_age %>%
  rownames_to_column(var = "LabID") %>%
  dplyr::filter(LabID %in% ringed_metadata_age$LabID) %>%
  column_to_rownames(var = "LabID")

# RERUN Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"


robust_aitch_dist_bearded_age <- vegdist(asv_mat_bearded_age, method = "robust.aitchison")
head(rownames(robust_aitch_dist_bearded_age))

bearded_permanova_aitch_age <- adonis2(robust_aitch_dist_bearded_age ~ Age_Group, 
                                         data = bearded_metadata_age, 
                                         permutations = 999)

print(bearded_permanova_aitch_age) # p-value = 0.733 not significant


robust_aitch_dist_ringed_age <- vegdist(asv_mat_ringed_age, method = "robust.aitchison")
head(rownames(robust_aitch_dist_ringed_age))

ringed_permanova_aitch_age <- adonis2(robust_aitch_dist_ringed_age ~ Age_Group, 
                                         data = ringed_metadata_age, 
                                         permutations = 999)

print(ringed_permanova_aitch_age) # p-value = 0.401 significant



# metadata_beluga_Locale <- data.frame(sample_data(ps.major)) %>%
#   dplyr::filter(LabID %in% labels(robust_aitch_dist_beluga_Locale)) %>%
#   dplyr::arrange(match(LabID, labels(robust_aitch_dist_beluga_Locale)))
# 
# stopifnot(
#   identical(metadata$LabID, labels(robust_aitch_dist))
# )

# none significant 

# ------------------------------------------------------------------------------
# PCoA ORDIANTION PLOTS
# ------------------------------------------------------------------------------

# Host palette: retains your established ggplot default host colors
host_colors <- c(
  "bearded seal" = "#F8766D",
  "beluga whale" = "#00BA38",
  "ringed seal" = "#619CFF"
)

# Season palette
season_colors <- c(
  "Winter" = "#4E79A7",
  "Spring" = "#59A14F",
  "Summer" = "#E15759",
  "Autumn" = "#B07AA1"
)

# Locale palette
locale_colors <- c(
  "Arctic" = "#4E79A7",
  "South Bering" = "#59A14F",
  "Cook Inlet" = "#E15759"
)

# Age-group palette
age_colors <- c(
  "Non-pup" = "#8C6D31",
  "Pup" = "#E377C2"
)

# ------------------------------------------------------------------------------
# FUNCTION: CREATE PCoA SCORES FROM A DISTANCE MATRIX
# ------------------------------------------------------------------------------

make_pcoa_df <- function(dist_obj, metadata_df, grouping_var, sample_id = "LabID") {
  
  # Make sample ID a plain character vector
  metadata_df <- metadata_df %>%
    dplyr::mutate(
      "{sample_id}" := as.character(.data[[sample_id]])
    ) %>%
    dplyr::filter(.data[[sample_id]] %in% labels(dist_obj)) %>%
    dplyr::arrange(match(.data[[sample_id]], labels(dist_obj)))
  
  # Ensure exact alignment between metadata and distance labels
  stopifnot(
    identical(
      as.character(metadata_df[[sample_id]]),
      labels(dist_obj)
    )
  )
  
  # PCoA / classical multidimensional scaling
  # add = TRUE applies a correction if the distance produces negative eigenvalues.
  pcoa_out <- cmdscale(
    dist_obj,
    k = 2,
    eig = TRUE,
    add = TRUE
  )
  
  # Calculate explained variation from positive eigenvalues only.
  eig_positive <- pcoa_out$eig[pcoa_out$eig > 0]
  
  axis1_pct <- round(
    100 * eig_positive[1] / sum(eig_positive),
    1
  )
  
  axis2_pct <- round(
    100 * eig_positive[2] / sum(eig_positive),
    1
  )
  
  # Construct an ordination data frame and join relevant metadata.
  pcoa_df <- as.data.frame(pcoa_out$points) %>%
    tibble::rownames_to_column(var = sample_id) %>%
    dplyr::rename(
      PCoA1 = V1,
      PCoA2 = V2
    ) %>%
    dplyr::left_join(
      metadata_df,
      by = sample_id
    )
  
  list(
    data = pcoa_df,
    axis1_pct = axis1_pct,
    axis2_pct = axis2_pct,
    grouping_var = grouping_var
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: CREATE A PCoA PLOT
# ------------------------------------------------------------------------------

plot_pcoa <- function(
    pcoa_list,
    grouping_var,
    colors,
    legend_title,
    title,
    subtitle = NULL,
    ellipse = TRUE
) {
  
  grouping_var <- rlang::ensym(grouping_var)
  
  p <- ggplot(
    pcoa_list$data,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = !!grouping_var
    )
  )
  
  # Ellipses require at least three samples per group. They are descriptive only.
  if (ellipse) {
    p <- p +
      stat_ellipse(
        aes(colour = !!grouping_var),
        type = "t",
        level = 0.95,
        linewidth = 0.7,
        show.legend = FALSE
      )
  }
  
  p +
    geom_point(
      shape = 21,
      colour = "black",
      stroke = 0.35,
      size = 3,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    scale_colour_manual(
      values = colors,
      drop = FALSE,
      na.value = "grey70"
    ) +
    labs(
      x = paste0("PCoA 1 (", pcoa_list$axis1_pct, "%)"),
      y = paste0("PCoA 2 (", pcoa_list$axis2_pct, "%)"),
      fill = legend_title,
      title = title,
      subtitle = subtitle
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}

# ==============================================================================
# PCoA BY HOST
# ==============================================================================

# metadata is already aligned to robust_aitch_dist after your stopifnot above.
# Re-align defensively in case you re-ran earlier chunks out of order.

metadata_host_pcoa <- metadata %>%
  dplyr::filter(LabID %in% labels(robust_aitch_dist)) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    Predator = factor(
      Predator,
      levels = c("bearded seal", "beluga whale", "ringed seal")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist)))

stopifnot(
  identical(
    metadata_host_pcoa$LabID,
    labels(robust_aitch_dist)
  )
)

pcoa_host <- make_pcoa_df(
  dist_obj = robust_aitch_dist,
  metadata_df = metadata_host_pcoa,
  grouping_var = "Predator"
)

p_pcoa_host <- plot_pcoa(
  pcoa_list = pcoa_host,
  grouping_var = Predator,
  colors = host_colors,
  legend_title = "Host",
  title = "Robust Aitchison PCoA by host",
  subtitle = "PERMANOVA: R² = 0.159, p = 0.001"
)

p_pcoa_host



# Optional centroid overlay for host plot
host_centroids <- pcoa_host$data %>%
  dplyr::group_by(Predator) %>%
  dplyr::summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  )

p_pcoa_host <- p_pcoa_host +
  geom_point(
    data = host_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = Predator
    ),
    inherit.aes = FALSE,
    shape = 23,
    colour = "black",
    stroke = 0.9,
    size = 5
  )

p_pcoa_host

# ==============================================================================
# PCoA BY SEASON WITHIN EACH HOST
# ==============================================================================

# ------------------------------------------------------------------------------
# BEARDED SEAL: SEASON
# ------------------------------------------------------------------------------

bearded_metadata_season <- bearded_metadata_season %>%
  dplyr::filter(
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    season = factor(
      season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_bearded_season)))

stopifnot(
  identical(
    bearded_metadata_season$LabID,
    labels(robust_aitch_dist_bearded_season)
  )
)

pcoa_bearded_season <- make_pcoa_df(
  dist_obj = robust_aitch_dist_bearded_season,
  metadata_df = bearded_metadata_season,
  grouping_var = "season"
)

p_pcoa_bearded_season <- plot_pcoa(
  pcoa_list = pcoa_bearded_season,
  grouping_var = season,
  colors = season_colors,
  legend_title = "Season",
  title = "Bearded seal robust Aitchison PCoA by season",
  subtitle = "PERMANOVA: R² = 0.081, p = 0.160",
  ellipse = TRUE
)

p_pcoa_bearded_season

# ------------------------------------------------------------------------------
# BELUGA WHALE: SEASON
# ------------------------------------------------------------------------------

beluga_metadata_season <- beluga_metadata_season %>%
  dplyr::filter(
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    season = factor(
      season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_beluga_season)))

stopifnot(
  identical(
    beluga_metadata_season$LabID,
    labels(robust_aitch_dist_beluga_season)
  )
)

pcoa_beluga_season <- make_pcoa_df(
  dist_obj = robust_aitch_dist_beluga_season,
  metadata_df = beluga_metadata_season,
  grouping_var = "season"
)

# Spring has only one beluga sample and winter has three samples; omit ellipses
# because 95% ellipse estimates are unstable or unavailable in such small groups.
p_pcoa_beluga_season <- plot_pcoa(
  pcoa_list = pcoa_beluga_season,
  grouping_var = season,
  colors = season_colors,
  legend_title = "Season",
  title = "Beluga robust Aitchison PCoA by season",
  subtitle = "PERMANOVA: R² = 0.088, p = 0.234",
  ellipse = FALSE
)

p_pcoa_beluga_season

# ------------------------------------------------------------------------------
# RINGED SEAL: SEASON
# ------------------------------------------------------------------------------

ringed_metadata_season <- ringed_metadata_season %>%
  dplyr::filter(
    !is.na(season),
    trimws(as.character(season)) != ""
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    season = factor(
      season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_ringed_season)))

stopifnot(
  identical(
    ringed_metadata_season$LabID,
    labels(robust_aitch_dist_ringed_season)
  )
)

pcoa_ringed_season <- make_pcoa_df(
  dist_obj = robust_aitch_dist_ringed_season,
  metadata_df = ringed_metadata_season,
  grouping_var = "season"
)

p_pcoa_ringed_season <- plot_pcoa(
  pcoa_list = pcoa_ringed_season,
  grouping_var = season,
  colors = season_colors,
  legend_title = "Season",
  title = "Ringed seal robust Aitchison PCoA by season",
  subtitle = "PERMANOVA: R² = 0.083, p = 0.120",
  ellipse = TRUE
)

p_pcoa_ringed_season

# Combine season PCoA plots
pcoa_plots_season <- (
  p_pcoa_bearded_season +
    p_pcoa_beluga_season +
    p_pcoa_ringed_season
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

pcoa_plots_season

ggsave(
  "Deliverables/ubiome/betadiv/AITCH_PCoA_by-season_each-host.png",
  plot = pcoa_plots_season,
  width = 18,
  height = 6,
  units = "in",
  dpi = 300
)

# ==============================================================================
# PCoA BY LOCALE WITHIN EACH HOST
# ==============================================================================

# Do not plot bearded seal locale effects:
# all bearded-seal samples are Arctic, so there is no within-host locale contrast.

# ------------------------------------------------------------------------------
# BELUGA WHALE: LOCALE
# ------------------------------------------------------------------------------

beluga_metadata_Locale <- beluga_metadata_Locale %>%
  dplyr::filter(
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    Locale = factor(
      Locale,
      levels = c("Arctic", "South Bering", "Cook Inlet")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_beluga_Locale)))

stopifnot(
  identical(
    beluga_metadata_Locale$LabID,
    labels(robust_aitch_dist_beluga_Locale)
  )
)

pcoa_beluga_Locale <- make_pcoa_df(
  dist_obj = robust_aitch_dist_beluga_Locale,
  metadata_df = beluga_metadata_Locale,
  grouping_var = "Locale"
)

# PERMDISP was significant for beluga locale (p = 0.0278), so this plot should
# be interpreted as showing both potential centroid separation and unequal spread.
p_pcoa_beluga_Locale <- plot_pcoa(
  pcoa_list = pcoa_beluga_Locale,
  grouping_var = Locale,
  colors = locale_colors,
  legend_title = "Locale",
  title = "Beluga robust Aitchison PCoA by locale",
  subtitle = "PERMANOVA: p = 0.006; PERMDISP: p = 0.028"
)

p_pcoa_beluga_Locale

# ------------------------------------------------------------------------------
# RINGED SEAL: LOCALE
# ------------------------------------------------------------------------------

ringed_metadata_Locale <- ringed_metadata_Locale %>%
  dplyr::filter(
    !is.na(Locale),
    trimws(as.character(Locale)) != ""
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    Locale = factor(
      Locale,
      levels = c("Arctic", "South Bering", "Cook Inlet")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_ringed_Locale)))

stopifnot(
  identical(
    ringed_metadata_Locale$LabID,
    labels(robust_aitch_dist_ringed_Locale)
  )
)

pcoa_ringed_Locale <- make_pcoa_df(
  dist_obj = robust_aitch_dist_ringed_Locale,
  metadata_df = ringed_metadata_Locale,
  grouping_var = "Locale"
)

# Ringed seals have only two South Bering samples; do not draw confidence
# ellipses because the ellipse is not stable or informative at n = 2.
p_pcoa_ringed_Locale <- plot_pcoa(
  pcoa_list = pcoa_ringed_Locale,
  grouping_var = Locale,
  colors = locale_colors,
  legend_title = "Locale",
  title = "Ringed seal robust Aitchison PCoA by locale",
  subtitle = "PERMANOVA: p = 0.911; South Bering n = 2",
  ellipse = FALSE
)

p_pcoa_ringed_Locale

# Combine locale PCoA plots
pcoa_plots_Locale <- (
  p_pcoa_beluga_Locale #+
    # p_pcoa_ringed_Locale
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

pcoa_plots_Locale

ggsave(
  "Deliverables/ubiome/betadiv/AITCH_PCoA_by-locale_each-host.png",
  plot = pcoa_plots_Locale,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)

# ==============================================================================
# PCoA BY AGE GROUP WITHIN EACH SEAL HOST
# ==============================================================================

# Metadata_age was already filtered to remove Unknown age records before the
# age-specific robust Aitchison distances were calculated.

# ------------------------------------------------------------------------------
# BEARDED SEAL: AGE GROUP
# ------------------------------------------------------------------------------

bearded_metadata_age <- bearded_metadata_age %>%
  dplyr::filter(
    !is.na(Age_Group),
    trimws(as.character(Age_Group)) != "",
    Age_Group %in% c("Non-pup", "Pup")
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    Age_Group = factor(
      Age_Group,
      levels = c("Non-pup", "Pup")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_bearded_age)))

stopifnot(
  identical(
    bearded_metadata_age$LabID,
    labels(robust_aitch_dist_bearded_age)
  )
)

pcoa_bearded_age <- make_pcoa_df(
  dist_obj = robust_aitch_dist_bearded_age,
  metadata_df = bearded_metadata_age,
  grouping_var = "Age_Group"
)

p_pcoa_bearded_age <- plot_pcoa(
  pcoa_list = pcoa_bearded_age,
  grouping_var = Age_Group,
  colors = age_colors,
  legend_title = "Age group",
  title = "Bearded seal robust Aitchison PCoA by age group",
  subtitle = "PERMANOVA: p = 0.733"
)

p_pcoa_bearded_age

# ------------------------------------------------------------------------------
# RINGED SEAL: AGE GROUP
# ------------------------------------------------------------------------------

ringed_metadata_age <- ringed_metadata_age %>%
  dplyr::filter(
    !is.na(Age_Group),
    trimws(as.character(Age_Group)) != "",
    Age_Group %in% c("Non-pup", "Pup")
  ) %>%
  dplyr::mutate(
    LabID = as.character(LabID),
    Age_Group = factor(
      Age_Group,
      levels = c("Non-pup", "Pup")
    )
  ) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_ringed_age)))

stopifnot(
  identical(
    ringed_metadata_age$LabID,
    labels(robust_aitch_dist_ringed_age)
  )
)

pcoa_ringed_age <- make_pcoa_df(
  dist_obj = robust_aitch_dist_ringed_age,
  metadata_df = ringed_metadata_age,
  grouping_var = "Age_Group"
)

p_pcoa_ringed_age <- plot_pcoa(
  pcoa_list = pcoa_ringed_age,
  grouping_var = Age_Group,
  colors = age_colors,
  legend_title = "Age group",
  title = "Ringed seal robust Aitchison PCoA by age group",
  subtitle = "PERMANOVA: p = 0.401"
)

p_pcoa_ringed_age

# Combine age-group PCoA plots
pcoa_plots_age <- (
  p_pcoa_bearded_age +
    p_pcoa_ringed_age
) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

pcoa_plots_age

ggsave(
  "Deliverables/ubiome/betadiv/AITCH_PCoA_by-age_each-host.png",
  plot = pcoa_plots_age,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)

# ==============================================================================
# FINAL COMBINED PCoA FIGURE
# ==============================================================================

all_pcoa_plots <- (
  p_pcoa_host /
    pcoa_plots_season #/
    # pcoa_plots_Locale /
    # pcoa_plots_age
) &
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

all_pcoa_plots

ggsave(
  "Deliverables/ubiome/betadiv/AITCH_ubiome-majorclass_PCoA.png",
  plot = all_pcoa_plots,
  width = 20,
  height = 24,
  units = "in",
  dpi = 300
)

# ==============================================================================
# END PCoA SECTION
# 