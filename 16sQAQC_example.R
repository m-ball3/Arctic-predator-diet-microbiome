library(tidyverse)
library(phyloseq)

### load phyloseq --------------------------------------------------------------

load("ps.16s.Rdata")
nsamples(ps.16s)
# replicates are removed already

## QAQC samples -----------------------

# Filter out anything not in Actinopteri
ps.16s_fish <- subset_taxa(ps.16s, Class == "Actinopteri")

nsamples(ps.16s_fish)

mean(sample_sums(ps.16s))
min(sample_sums(ps.16s))
median(sample_sums(ps.16s))
hist(sample_sums(ps.16s), breaks = seq(0,140000,500))

# Remove samples with total abundance < 100
ps.16s_filt <- prune_samples(sample_sums(ps.16s_fish) >= 100, ps.16s)

nsamples(ps.16s_filt)

mean(sample_sums(ps.16s_filt))
min(sample_sums(ps.16s_filt))
median(sample_sums(ps.16s_filt))
hist(sample_sums(ps.16s_filt), breaks = seq(0,140000,500))

## QAQC taxa -------------------------

## Merge by species
ps.16s.sp <- tax_glom(ps.16s, "Species.y", NArm = FALSE)

# Convert to proportion
ps16s_rel <- transform_sample_counts(ps.16s.sp, 
                                     function(x) {x / sum(x)})


#must be at least 1% of diet in 1 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps16s_rel, f1, A=1)

ps.sp_minor <- prune_taxa(lowcount.filt, ps16s_rel)

#must be at least 1% of diet in 4 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps16s_rel, f1, A=4)

ps.sp_major <- prune_taxa(lowcount.filt, ps16s_rel)
