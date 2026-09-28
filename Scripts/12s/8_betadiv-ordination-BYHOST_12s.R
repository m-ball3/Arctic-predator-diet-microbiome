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
load("./Scripts/12s/rdata/ps.12s.merged.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.12s.merged), "matrix")
head(rownames(asv_mat))
# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------
# formatting for metadata
metadata <- data.frame(sample_data(ps.12s.merged))

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

print(permanova_aitch_host) # p-value = 0.04 

metadata <- data.frame(sample_data(ps.12s.merged)) %>%
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
# 
# # PERMDIST
# ## this box plot will show within-host beta diversity, 
# ## or dispersion: whether samples within each predator group are similarly variable.
# 
# # creates a usable column from robust aitchison distances
# ## disp$distances is the distance from each sample to its host-group centroid
# 
# disp_plot_df <- metadata %>%
#   dplyr::mutate(
#     Distance_to_centroid = disp$distances
#   )
# 
# stopifnot(
#   nrow(disp_plot_df) == length(disp$distances),
#   identical(as.character(disp_plot_df$Predator), as.character(disp$group))
# )
# # Larger values mean samples are more dispersed/heterogeneous within that host group
# box_plot_host <- ggplot(
#   disp_plot_df,
#   aes(
#     x = Predator,
#     y = Distance_to_centroid,
#     fill = Predator
#   )
# ) +
#   geom_boxplot(
#     outlier.shape = NA,
#     alpha = 0.7,
#     colour = "black"
#   ) +
#   geom_jitter(
#     width = 0.15,
#     shape = 21,
#     colour = "black",
#     size = 2,
#     alpha = 0.8,
#     stroke = 0.35
#   ) +
#   labs(
#     x = "Host",
#     y = "Distance to group centroid",
#     title = "Multivariate dispersion by host"
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none"
#   )
# 
# box_plot_host
# 
# # WITHIN VS BETWEEN-HOST DISTANCES
# # creates a longform df of robust aitchison distances
# pairwise_dist_df <- as.matrix(robust_aitch_dist) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("LabID_1") %>%
#   tidyr::pivot_longer(
#     cols = -LabID_1,
#     names_to = "LabID_2",
#     values_to = "Robust_Aitchison"
#   ) %>%
#   dplyr::filter(LabID_1 < LabID_2) %>%
#   dplyr::left_join(
#     metadata %>%
#       dplyr::select(LabID, Predator) %>%
#       dplyr::rename(
#         LabID_1 = LabID,
#         Predator_1 = Predator
#       ),
#     by = "LabID_1"
#   ) %>%
#   dplyr::left_join(
#     metadata %>%
#       dplyr::select(LabID, Predator) %>%
#       dplyr::rename(
#         LabID_2 = LabID,
#         Predator_2 = Predator
#       ),
#     by = "LabID_2"
#   ) %>%
#   dplyr::mutate(
#     comparison = dplyr::case_when(
#       Predator_1 == Predator_2 ~ paste("Within", Predator_1),
#       TRUE ~ paste(
#         pmin(Predator_1, Predator_2),
#         "vs",
#         pmax(Predator_1, Predator_2)
#       )
#     )
#   )
# 
# # plots
# pairwise_dist_df %>%
#   ggplot(
#     aes(
#       x = comparison,
#       y = Robust_Aitchison,
#       fill = comparison
#     )
#   ) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.75) +
#   geom_jitter(
#     width = 0.15,
#     alpha = 0.15,
#     size = 1
#   ) +
#   labs(
#     x = NULL,
#     y = "Pairwise robust Aitchison distance",
#     title = "Within- and between-host microbiome dissimilarity"
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(
#       angle = 35,
#       hjust = 1
#     ),
#     legend.position = "none"
#   )
# 
# # some are significantly different. 
# # bearded seals and beluga whales (BH adjusted p-value = 0.00015)
# # ringed seals and beluga whales (BH adjusted p-value = 0.00015)
# # bearded seals and ringed seals are not significantly different (BH adjusted p-value = 0.0519)

# ------------------------------------------------------------------------------
# PCoA ORDIANTION PLOTS
# ------------------------------------------------------------------------------

# Host palette: retains your established ggplot default host colors
host_colors <- c(
  "bearded seal" = "#F8766D",
  "beluga whale" = "#00BA38",
  "ringed seal" = "#619CFF"
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
  subtitle = "PERMANOVA: R² = 0.0503, p = 0.04"
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

ggsave(
  "Deliverables/12s/betadiv/AITCH_12s_PCoA.png",
  plot = p_pcoa_host,
  width = 20,
  height = 24,
  units = "in",
  dpi = 300
)

# ==============================================================================
# END PCoA SECTION
# 