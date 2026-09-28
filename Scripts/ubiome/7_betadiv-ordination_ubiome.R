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
load("./Scripts/ubiome/rdata/ubiome.ps-major.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.major), "matrix")
head(rownames(asv_mat))
# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------
# formatting for metadata
metadata <- data.frame(sample_data(ps.major))

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

metadata <- data.frame(sample_data(ps.major)) %>%
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

# Samples in metadata but not in asv mat
setdiff(rownames(metadata_season), rownames(asv_mat))

# Samples in asv mat but not in metadata
setdiff(rownames(asv_mat), rownames(metadata_season))

# cleans avs table based on metadata_seasons
asv_mat_season <- as.data.frame(asv_mat) %>% 
  rownames_to_column(var = "row_id") %>% 
  dplyr::filter(!row_id %in% c("WADE-003-013", "WADE-003-148", "WADE-003-149")) %>% 
  column_to_rownames(var = "row_id")

# RERUN Robust Aitchison distance
## Robust method only log-ratio tranforms non-zeros
## thus, avoiding the need to create data + 1 
### therefore, it avoids biasing by creating "pseudocounts"
robust_aitch_dist_season <- vegdist(asv_mat_season, method = "robust.aitchison")
head(rownames(robust_aitch_dist_season))

permanova_aitch_season <- adonis2(robust_aitch_dist_season ~ season, 
                           data = metadata_season, 
                           permutations = 999)

print(permanova_aitch_season) # p-value = 0.097 not significant

metadata_season <- data.frame(sample_data(ps.major)) %>%
  dplyr::filter(LabID %in% labels(robust_aitch_dist_season)) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist_season)))

stopifnot(
  identical(metadata_season$LabID, labels(robust_aitch_dist_season))
)

set.seed(123)

pairwise_aitch_season <- pairwise_permanova(
  dist_obj = robust_aitch_dist_season,
  metadata = metadata_season,
  group_var = "season",
  sample_id = "LabID",
  permutations = 9999,
  p_adjust = "BH"
)

pairwise_aitch_season

# tests for dispersion
disp <- vegan::betadisper(
  robust_aitch_dist_season,
  group = metadata_season$season
)

set.seed(123)
vegan::permutest(disp, permutations = 9999) # p-value = 0.2476
## Multivariate dispersion did not differ significantly among host groups


# No significant differences between seasons


# BY LOCATION--------------------------------------------------------------------------------
permanova_aitch_location <- adonis2(robust_aitch_dist ~ Location, 
                                data = metadata, 
                                permutations = 999)

print(permanova_aitch_location) # p-value = 0.001

metadata <- data.frame(sample_data(ps.major)) %>%
  dplyr::filter(LabID %in% labels(robust_aitch_dist)) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist)))

stopifnot(
  identical(metadata$LabID, labels(robust_aitch_dist))
)


set.seed(123)

pairwise_aitch_location <- pairwise_permanova(
  dist_obj = robust_aitch_dist,
  metadata = metadata,
  group_var = "Location",
  sample_id = "LabID",
  permutations = 9999,
  p_adjust = "BH"
)

location <- as.data.frame(pairwise_aitch_location)

metadata %>%
  dplyr::count(Location, sort = TRUE)

## NOT ENOUGH SAMPLES IN EACH LOCATION TO CALCULATE THIS
# # tests for dispersion
# disp <- vegan::betadisper(
#   robust_aitch_dist,
#   group = metadata$Location
# )
# 
# set.seed(123)
# vegan::permutest(disp, permutations = 9999) 


# some locations were significantly different. Check:
# View(location)

# BY LOCALE----------------------------------------------------------------------------
permanova_aitch_locale <- adonis2(robust_aitch_dist ~ Locale, 
                                    data = metadata, 
                                    permutations = 999)

print(permanova_aitch_locale) # p-value = 0.001

metadata <- data.frame(sample_data(ps.major)) %>%
  dplyr::filter(LabID %in% labels(robust_aitch_dist)) %>%
  dplyr::arrange(match(LabID, labels(robust_aitch_dist)))

stopifnot(
  identical(metadata$LabID, labels(robust_aitch_dist))
)


set.seed(123)

pairwise_aitch_locale <- pairwise_permanova(
  dist_obj = robust_aitch_dist,
  metadata = metadata,
  group_var = "Locale",
  sample_id = "LabID",
  permutations = 9999,
  p_adjust = "BH"
)

pairwise_aitch_locale

# tests for dispersion
disp <- vegan::betadisper(
  robust_aitch_dist,
  group = metadata$Locale
)

set.seed(123)
vegan::permutest(disp, permutations = 9999) # p-value = 0.4162
## Multivariate dispersion did not differ significantly among locale groups
## Therefore, the significant PERMANOVA result is not accompanied by formal
## evidence of heterogeneous dispersion at alpha = 0.05.


# all three locales was significantly different
## Arctic was significantly different than South Bering (p = 0.0003)
## Cook inlet was significantly different than South Bering (p = 0.0275)
## Arctic was significantly different than Cook inlet (p = 0.0275)

# BY AGE------------------------------------------------------------------------------
permanova_aitch_age <- adonis2(robust_aitch_dist_age ~ Age_Group, 
                                  data = metadata_age, 
                                  permutations = 999)

print(permanova_aitch_age) # p-value = 0.063

set.seed(123)

pairwise_aitch_age <- pairwise_permanova(
  dist_obj = robust_aitch_dist_age,
  metadata = metadata_age,
  group_var = "Age_Group",
  sample_id = "LabID",
  permutations = 9999,
  p_adjust = "BH"
)

pairwise_aitch_age

# OUT O FBOUNDS
# # tests for dispersion
# disp <- vegan::betadisper(
#   robust_aitch_dist,
#   group = metadata$Age_Group
# )
# 
# set.seed(123)
# vegan::permutest(disp, permutations = 9999) # p-value = 0.4162


# no significant differences between age classes


# # OPTIONAL PCA -------------------------------------------------------------------------------------------

# rCLR transform (manually creates the matrix that vegdist creates)
otu_rclr <- decostand(asv_mat, method = "rclr")
otu_rclr_age <- decostand(asv_mat_age, method = "rclr")

# Runs PCA on rCLR-transformed data
pca_res <- prcomp(otu_rclr, center = TRUE, scale. = FALSE)
pca_res

pca_res_age <- prcomp(otu_rclr_age, center = TRUE, scale. = FALSE)
pca_res_age

# Scores (sample coordinates in PCA space)
pca_scores <- pca_res$x
pca_scores

pca_scores_age <- pca_res_age$x
pca_scores_age

# Loadings (taxon contributions)
pca_loadings <- pca_res$rotation
pca_loadings

pca_loadings_age <- pca_res_age$rotation
pca_loadings_age


# PCA scores (samples × PCs)
pca_scores <- as.data.frame(pca_res$x)
pca_scores_age <- as.data.frame(pca_res_age$x)

# Add sample IDs as a column
pca_scores$LabID <- rownames(pca_scores)
pca_scores_age$LabID <- rownames(pca_scores_age)

# Join scores with metadata
pca_plot_df <- pca_scores %>%
  dplyr::left_join(metadata, by = "LabID")

var_explained <- summary(pca_res)$importance[2, ]  # proportion per PC [web:244]
pc1_pct <- round(var_explained[1] * 100, 1)
pc2_pct <- round(var_explained[2] * 100, 1)

pca_plot_df_age <- pca_scores_age %>%
  dplyr::left_join(metadata_age, by = "LabID")

var_explained_age <- summary(pca_res_age)$importance[2, ]  # proportion per PC [web:244]
pc1_pct_age <- round(var_explained_age[1] * 100, 1)
pc2_pct_age <- round(var_explained_age[2] * 100, 1)

# BY HOST
p_pca_host <- ggplot(pca_plot_df,
                        aes(x = PC1, y = PC2,
                            color = Predator)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95)+
  # facet_wrap(~ Sample_index, ncol = 4) +
  labs(
    x = paste0("PC1 (", pc1_pct, "% variance)"),
    y = paste0("PC2 (", pc2_pct, "% variance)"),
    color = "Host",
    title = "Aitchison Distances by host"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

p_pca_host

### Biological interpretation
# Beta-diversity results support the conclusion that bearded-seal, ringed-seal, 
# and beluga-whale samples have distinct overall microbiome compositions under 
# the robust Aitchison framework. The PCA visualization is consistent with that 
# result, but it also shows substantial within-host variation and overlap—so 
# host species explains an important, but not complete, component of microbiome variation.

# BY SEASON
p_pca_season <- ggplot(pca_plot_df,
                     aes(x = PC1, y = PC2,
                         color = season)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95)+
  # facet_wrap(~ Sample_index, ncol = 4) +
  labs(
    x = paste0("PC1 (", pc1_pct, "% variance)"),
    y = paste0("PC2 (", pc2_pct, "% variance)"),
    color = "season", 
    title = "Aitchison Distances by season"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

p_pca_season

### Biological interpretation

# BY LOCATION
p_pca_location <- ggplot(pca_plot_df,
                     aes(x = PC1, y = PC2,
                         color = Location)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95)+
  # facet_wrap(~ Sample_index, ncol = 4) +
  labs(
    x = paste0("PC1 (", pc1_pct, "% variance)"),
    y = paste0("PC2 (", pc2_pct, "% variance)"),
    color = "location", 
    title = "Aitchison Distances by location"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

p_pca_location

### Biological interpretation

# BY LOCALE
p_pca_locale <- ggplot(pca_plot_df,
                     aes(x = PC1, y = PC2,
                         color = Locale)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95)+
  # facet_wrap(~ Sample_index, ncol = 4) +
  labs(
    x = paste0("PC1 (", pc1_pct, "% variance)"),
    y = paste0("PC2 (", pc2_pct, "% variance)"),
    color = "locale",
    title = "Aitchison Distances by locale"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

p_pca_locale

### Biological interpretation

# BY AGE
p_pca_age <- ggplot(pca_plot_df_age,
                     aes(x = PC1, y = PC2,
                         color = Age_Group)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(level = 0.95)+
  # facet_wrap(~ Sample_index, ncol = 4) +
  labs(
    x = paste0("PC1 (", pc1_pct, "% variance)"),
    y = paste0("PC2 (", pc2_pct, "% variance)"),
    color = "age group", 
    title = "Aitchison Distances by age group"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 12),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

p_pca_age

### Biological interpretation


all_plots <- (p_pca_host + p_pca_season) / (p_pca_location + p_pca_locale) / p_pca_age &theme(
  text = element_text(size = 30),
  axis.title = element_text(size = 30),
  axis.text = element_text(size = 30),
  plot.title = element_text(size = 30),
  legend.title = element_text(size = 30),
  legend.text = element_text(size = 30)
)

ggsave("Deliverables/ubiome/betadiv/AITCH_ubiome-majorclass.png", plot = all_plots, width = 40, height = 35, units = "in", dpi = 300)
