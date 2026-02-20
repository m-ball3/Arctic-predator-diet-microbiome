library(tidyverse)
library(phyloseq)

### load phyloseq --------------------------------------------------------------

load("ps.12s.Rdata")
nsamples(ps.12s)
# replicates are removed already

hist(sample_sums(ps.12s), breaks = seq(0,2000000,500))

## QAQC samples -----------------------

# Filter out anything not in Actinopteri
ps.12s_fish <- subset_taxa(ps.12s, Class == "Actinopteri")

nsamples(ps.12s_fish)

mean(sample_sums(ps.12s_fish))
min(sample_sums(ps.12s_fish))
median(sample_sums(ps.12s_fish))
hist(sample_sums(ps.12s_fish), breaks = seq(0,2000000,100))

# Remove samples with total abundance < 100
ps.12s_filt <- prune_samples(sample_sums(ps.12s_fish) >= 100, ps.12s)

nsamples(ps.12s_filt)

mean(sample_sums(ps.12s_filt))
min(sample_sums(ps.12s_filt))
median(sample_sums(ps.12s_filt))
hist(sample_sums(ps.12s_filt), breaks = seq(0,2000000,500))

## QAQC taxa -------------------------

## Merge by species
ps.12s.sp <- tax_glom(ps.12s_filt, "Species.y", NArm = FALSE)

# Convert to proportion
ps12s_rel <- transform_sample_counts(ps.12s.sp, 
                                     function(x) {x / sum(x)}) %>% 
  prune_samples(sample_sums(ps.12s.sp) >= 0, .)

plot_bar(ps12s_rel, x = "LabID", fill = "Species.y") + 
  facet_wrap(~Predator, ncol = 1, scales = "free_x")

#must be at least 1% of diet in 1 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps12s_rel, f1, A=1)

ps12s.sp_minor <- prune_taxa(lowcount.filt, ps12s_rel)

plot_bar(ps12s.sp_minor, x = "LabID", fill = "Species.y") + 
  facet_wrap(~Predator, ncol = 1, scales = "free_x")

#must be at least 1% of diet in 4 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps12s_rel, f1, A=4)

ps12s.sp_major <- prune_taxa(lowcount.filt, ps12s_rel) #%>% 
  #prune_samples(sample_sums(.) > 0, .)

plot_bar(ps12s.sp_major, x = "LabID", fill = "Species.y") + 
  facet_wrap(~Predator, ncol = 1, scales = "free_x")
