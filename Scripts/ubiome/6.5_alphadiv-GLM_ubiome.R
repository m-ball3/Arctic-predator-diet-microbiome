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
library(MASS)
library(DHARMa)
library(emmeans)


getwd()
load("./Scripts/ubiome/rdata/ubiome.ps-major.Rdata")

# gets the otu table
asv_mat <- as(otu_table(ps.major), "matrix")

# ------------------------------------------------------------------------------
# Formats data for downstream analysis
# ------------------------------------------------------------------------------

# formatting for metadata
metadata <- data.frame(sample_data(ps.major))

# sets up Predator as a factor
metadata <- metadata %>%
  dplyr::mutate(Predator = factor(Predator))

class(metadata)

# optionally removes fin whale, since there is only one sample (two sub-samples of the same sample)
metadata <- metadata %>%
  dplyr::filter(Predator != "fin whale")

asv_mat <- asv_mat[
  !rownames(asv_mat) %in% c("WADE-003-150-A", "WADE-003-150-B", "neg1","mock1-782D","mock2-783D","neg2","neg2-2"),,
  drop = FALSE
]

# sanity checks
# samples_in_asv_not_metadata <- setdiff(
#   rownames(asv_mat),
#   metadata$LabID
# )
# 
# samples_in_asv_not_metadata

# Puts metadata into precisely the same order as the ASV table
metadata <- metadata %>%
  dplyr::filter(LabID %in% rownames(asv_mat)) %>%
  dplyr::slice(match(rownames(asv_mat), LabID))

# sanity check: all should return TRUE
identical(rownames(asv_mat), metadata$LabID)
nrow(asv_mat) == nrow(metadata)

# ------------------------------------------------------------------------------
# Adds diet data
# ------------------------------------------------------------------------------
load("./Scripts/12s/rdata/ps.12s.dietdat.Rdata")

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

# filters for only samples that have diet data
metadata <- metadata %>%
  dplyr::filter(!is.na(diet_observed))

# sees how many samples are left per host
metadata %>%
  dplyr::count(Predator)
# bearded seal 20
# beluga whale 26
# ringed seal 23

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

# ------------------------------------------------------------------------------
# Richness
# ------------------------------------------------------------------------------

# Robbins
## the ratio of singletons to total taxa


# Observed 
## the number of ASVs present
# ------------------------------------------------------------------------------
# Tests for normality
# ------------------------------------------------------------------------------

# Visual and statistical tests for normality

# Histogram of the relative abundance of chinook
# ?hist()
hist(metadata$Observed)


# Shapiro Wilkes test for normality
# makes Observed into a numeric vector
y <- metadata$Observed
range(y)

# Runs Shapiro-Wilk normality test
## Null hypothesis: the data are normal
## Alt hypothesis: the data are not normal
## alpha = 0.05

shapiro.test(y) #p = 0.2104 -> insufficient evidence to reject normality



# test homogeneity of variance FOR HOST
?leveneTest()
var <- leveneTest(Observed ~ Predator, metadata)
var # p-value = 0.4171 -> no evidence of unequal raw observed variance among predator groups

# Meets assumptions for a parametric model
## normal distrudfvdfvdei
# ------------------------------------------------------------------------------
# Fits simple models to test whether to use a Poisson or Negative Binomial
# ------------------------------------------------------------------------------

p.mod1  <- glm(Observed ~ Predator,
               family = poisson(link = "log"),
               data = metadata)

nb.mod1 <- glm.nb(Observed ~ Predator,
                  data = metadata)

# 1. Compare information criteria
AIC(p.mod1, nb.mod1) # essentially the same
BIC(p.mod1, nb.mod1) # essentially the same

# 2. Compare coefficient estimates / uncertainty
summary(p.mod1)
summary(nb.mod1)

# 3. Simulation-based residual diagnostics
set.seed(123)
p_res  <- simulateResiduals(p.mod1,  n = 1000)
nb_res <- simulateResiduals(nb.mod1, n = 1000)

plot(p_res)
plot(nb_res)

testDispersion(p_res)
testDispersion(nb_res)

testZeroInflation(p_res)
testZeroInflation(nb_res)


## negative binomial is a more complicated model and it does not add anything. '
# Moving forward with the Poisson

# ------------------------------------------------------------------------------
# Fits and diagnoses an additive vs interactive poisson model
# ------------------------------------------------------------------------------

# Poisson generalized linear model (additive)
p.mod2 <- glm(
  Observed ~ Predator + diet_observed,
  family = poisson(link = "log"),
  data = metadata
)
plot(p.mod2)

# Poisson-model diagnostics
set.seed(123)
p_res2 <- DHARMa::simulateResiduals(
  fittedModel = p.mod2,
  n = 1000
)

plot(p_res2)

DHARMa::testUniformity(p_res2)
DHARMa::testDispersion(p_res2)

# Poisson generalized linear model (interaction term)
p.mod3 <- glm(
  Observed ~ Predator*diet_observed,
  family = poisson(link = "log"),
  data = metadata
)
plot(p.mod3)

# Poisson-model diagnostics
set.seed(123)
p_res3 <- DHARMa::simulateResiduals(
  fittedModel = p.mod3,
  n = 1000
)

plot(p_res3)

DHARMa::testUniformity(p_res3)
DHARMa::testDispersion(p_res3)

# 1. Compare information criteria
AIC(p.mod2, p.mod3) # essentially the same
BIC(p.mod2, p.mod3) # essentially the same

# 2. Compare coefficient estimates / uncertainty
summary(p.mod2)
summary(p.mod3)

# 3. Anova
anova(p.mod2, p.mod3, test = "Chisq") # p = 0.026

# conclusion: the interaction model does contribute information, so I am moving 
# forward with it
## the interaction term allows for separate slopes per predator, which makes sense
## to me biologically (evidence for significant differences by host in other tests)

# ------------------------------------------------------------------------------
# Should read depth be entered as a random effect in this model?
## also, should read depth for the diet also be entered as a random effect??
# ------------------------------------------------------------------------------

# tests for correlations between library depth and observed richness
cor.test(
  metadata$Observed,
  metadata$n_reads,
  method = "spearman",
  exact = FALSE
)

# plots
ggplot(metadata, aes(x = n_reads, y = Observed, colour = Predator)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    colour = "black",
    linewidth = 0.8
  ) +
  scale_x_log10() +
  labs(
    x = "16S library size (reads; log10 scale)",
    y = "Observed ASV richness",
    colour = "Predator"
  ) +
  theme_classic()

# p-value = 0.03582; there is a negative trend (rho = -0.2531736), there is evidence
# that the read depth is correlated with the observed richness
# when looking at the plot, there are some possible outlier with high leverage
### interpretation = account for this in the model

# Tests for differences in read depth by host
metadata %>%
  dplyr::group_by(Predator) %>%
  dplyr::summarise(
    n = dplyr::n(),
    spearman_rho = cor(
      Observed, n_reads,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

# group specific tests
metadata %>%
  dplyr::group_by(Predator) %>%
  dplyr::group_modify(~ {
    out <- cor.test(
      .x$Observed,
      .x$n_reads,
      method = "spearman",
      exact = FALSE
    )
    
    tibble::tibble(
      n = nrow(.x),
      rho = unname(out$estimate),
      p_value = out$p.value
    )
  })



# there is evidence that the read depth by host is correlated with the observed richness

# ------------------------------------------------------------------------------
# GLMM (adds random effects for read depth to the poisson model with diet as an interaction term)
# ------------------------------------------------------------------------------

# Poisson generalized linear model (interaction term)
p.mod4 <- glm(
  Observed ~ Predator*diet_observed,
  family = poisson(link = "log"),
  data = metadata
)
plot(p.mod4)

# Poisson-model diagnostics
set.seed(123)
p_res3 <- DHARMa::simulateResiduals(
  fittedModel = p.mod3,
  n = 1000
)

plot(p_res3)

DHARMa::testUniformity(p_res3)
DHARMa::testDispersion(p_res3)

# 1. Compare information criteria
AIC(p.mod2, p.mod3) # essentially the same
BIC(p.mod2, p.mod3) # essentially the same

# 2. Compare coefficient estimates / uncertainty
summary(p.mod2)
summary(p.mod3)

# 3. Anova
anova(p.mod2, p.mod3, test = "Chisq") # p = 0.026

# conclusion: the additive model does contribute information, so I am moving 
# forward with it
## the interaction term allows for separate slopes per predator, which makes sense
## to me biologically (evidence for significant differences by host in other tests)


# ------------------------------------------------------------------------------
# 2. Mixed Model Analysis  FROM SWEENEY ET AL. 2023

# Mixed models here will all be carried out in MCMCglmm, 
# but your favourite package/ program should work just fine with some adaptation as well. 
# Feel free to contact me about syntax switches between algorithms. 

## 2A. Poisson GLMMs 

# * These models produce the main results of the manuscript and 
# use the raw read count per ASV per sample as the response with Poisson error families 
# * Below we also present alternative approaches for computing speed and adjustment 
# of questions asked by the model 

### Model set-up 

# * read abundance as response 
# * fixed effects should be relevant host factors that may affect sample variation 


# response: raw read abundance per ASV per sample
# fixed effects: 
## Predator = average community-wide difference in ASV abundance among hosts.
## Diet = 

# random effects: sample identity, ASV identity, and ASV-by-predator variation
# Random effects:
#   (1 | LabID)
#     Sample-specific baseline abundance across ASVs. This absorbs global
#     sample-to-sample count variation, including library size, but may also
#     reflect biological variation in total bacterial biomass/community load.
#
#   (1 | ASV)
#     Baseline abundance differences among ASVs. Large variance is expected:
#     some ASVs are consistently much more abundant than others.
#
#   (1 | ASV:Predator) ### DO WE NEED THIS INTERACTION? KEEP PRED FOR SURE
#     Taxon-specific variation among hosts. A non-zero variance indicates that
#     predator-associated differences are not uniform across ASVs: individual
#     ASVs can increase, decrease, or remain unchanged across hosts even if
#     the average fixed Predator effect is weak or non-significant.
#
# Distribution:
#   Negative binomial (nbinom2) accommodates overdispersion in raw ASV counts.

# ------------------------------------------------------------------------------
# Build one row per Sample × ASV.
# ------------------------------------------------------------------------------

# pivots asv_mat from wide form to longform
asv_long <- as.data.frame(asv_mat) %>%
  tibble::rownames_to_column("LabID") %>%
  tidyr::pivot_longer(
    cols = -LabID,
    names_to = "ASV",
    values_to = "abundance"
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(
        LabID,
        Predator,
        diet_observed,
        season
      ),
    by = "LabID"
  ) %>%
  dplyr::filter(
    !is.na(Predator),
    !is.na(diet_observed),
    !is.na(season)
  ) %>%
  dplyr::mutate(
    LabID = factor(LabID),
    ASV = factor(ASV),
    Predator = factor(Predator),
    season = factor(season),
    
    # Center/scale continuous predictors before interaction terms.
    diet_observed_z = as.numeric(scale(diet_observed))
  )

# sanity checks
stopifnot(all(asv_long$abundance >= 0))
stopifnot(all(asv_long$abundance == as.integer(asv_long$abundance)))

asv_long %>%
  dplyr::summarise(
    samples = dplyr::n_distinct(LabID),
    ASVs = dplyr::n_distinct(ASV),
    rows = dplyr::n(),
    zero_fraction = mean(abundance == 0)
  )

# fits a model
m_asv_predator <- glmmTMB::glmmTMB(
  abundance ~ Predator +
    (1 | LabID) +
    (1 | ASV),
   # (1 | ASV:Predator),
  family = glmmTMB::nbinom2,
  data = asv_long
)

m_asv_predator$sdr$pdHess
glmmTMB::diagnose(m_asv_predator)
summary(m_asv_predator)


# non significant predator effects signal that there is no community-wide direction by predator
## this does not mean that there is not differentially abundant taxa
### in fact, the variances for ASV and ASV:Predator indicate differentially abundant ASVs

# Therefore, significant host differences in richness may be driven by the
# presence or absence of low-abundance taxa, or by taxon-specific abundance 
# changes that do not occur uniformly across the community.


# ------------------------------------------------------------------------------
# ADDS DIET
## Does the number of distinct prey taxa detected in an animal’s diet relate 
## to its bacterial ASV abundance pattern?
# ------------------------------------------------------------------------------

# Models diet and host (additive)
m_asv_predator_diet1 <- glmmTMB::glmmTMB(
  abundance ~ diet_observed_z +
    (1 | LabID) +
    (1 | ASV) +
    (1 | ASV:Predator),
  family = glmmTMB::nbinom2,
  data = asv_long
)

m_asv_predator_diet1$sdr$pdHess
glmmTMB::diagnose(m_asv_predator_diet1)
summary(m_asv_predator_diet1)

# Models diet and host (interaction term)
m_asv_predator_diet2 <- glmmTMB::glmmTMB(
  abundance ~ Predator * diet_observed_z +
    (1 | LabID) +
    (1 | ASV) +
    (1 | ASV:Predator),
  family = glmmTMB::nbinom2,
  data = asv_long
)

m_asv_predator_diet2$sdr$pdHess
glmmTMB::diagnose(m_asv_predator_diet2)
summary(m_asv_predator_diet2)

anova(
  m_asv_predator_diet1,
  m_asv_predator_diet2
)

AIC(
  m_asv_predator_diet1,
  m_asv_predator_diet2
)

## evidence for the model with the interaction term, moving forward with that one. 

# the number of distinct prey taxa detected in an animal's diet does relate to the 
# bacterial ASV pattern, and this is stronger for belugas vs bearded seals. 



# NEED TO ADD A TEST FOR WHETHER LIBRARY SIZE AFFECTS DIET RICHNESS TO CAPTURE THAT
# EFFECT IN THE ABOVE MODEL!!!: 

metadata <- metadata %>%
  dplyr::left_join(
    metadata_12s %>%
      dplyr::select(
        LabID,
        diet_observed,
        diet_n_reads = n_reads
      ),
    by = "LabID"
  )

cor.test(
  metadata$diet_observed.y,
  metadata$diet_n_reads,
  method = "spearman",
  exact = FALSE
)

# p = 0.004, there is evidence that read depth matters. rho = .33
ggplot(
  metadata,
  aes(
    x = diet_n_reads,
    y = diet_observed.y,
    colour = Predator
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  scale_x_log10() +
  labs(
    x = "12S diet-library size (reads; log10 scale)",
    y = "Observed diet ASV richness",
    colour = "Predator"
  ) +
  theme_classic()

# PLOTS ------------------------------------------------------------------------

host_colors <- c(
  "bearded seal" = "#D55E00",
  "beluga whale" = "#009E73",
  "ringed seal"  = "#0072B2"
)

pred_grid <- tidyr::expand_grid(
  Predator = factor(
    levels(asv_long$Predator),
    levels = levels(asv_long$Predator)
  ),
  diet_observed_z = seq(
    min(asv_long$diet_observed_z, na.rm = TRUE),
    max(asv_long$diet_observed_z, na.rm = TRUE),
    length.out = 100
  )
)

# Population-level predictions:
# re.form = NA sets all random effects to zero.
pred_grid$predicted_abundance <- predict(
  m_asv_predator_diet2,
  newdata = pred_grid,
  type = "response",
  re.form = NA
)

# Plot
ggplot(
  pred_grid,
  aes(
    x = diet_observed_z,
    y = predicted_abundance,
    colour = Predator
  )
) +
  geom_line(
    linewidth = 1.4
  ) +
  scale_colour_manual(
    values = host_colors
  ) +
  labs(
    x = "Diet ASV richness (standardized)",
    y = "Predicted mean ASV read count",
    colour = "Predator"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )


# strong positive correlation. Add 12s reads as a covariate
asv_long <- asv_long %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(
        LabID,
        diet_n_reads
      ),
    by = "LabID"
  ) %>%
  dplyr::mutate(
    diet_log_reads_z = as.numeric(
      scale(log(diet_n_reads))
    )
  )

# Models the predator and diet interaction while controlling for 12s read depth
m_asv_predator_diet_depth <- glmmTMB::glmmTMB(
  abundance ~ Predator * diet_observed_z +
    diet_log_reads_z +
    (1 | LabID) +
    (1 | ASV) +
    (1 | ASV:Predator),
  family = glmmTMB::nbinom2,
  data = asv_long
)

m_asv_predator_diet_depth$sdr$pdHess
glmmTMB::diagnose(m_asv_predator_diet_depth)
summary(m_asv_predator_diet_depth)

anova(
  m_asv_predator_diet2,
  m_asv_predator_diet_depth
)

AIC(
  m_asv_predator_diet2,
  m_asv_predator_diet_depth
)

## adding the control for read depth did not really change the results from the model. 
### retain the model without the read depth covariate.

# Final ASV model: m_asv_predator_diet2
m_asv_final <- m_asv_predator_diet2
  