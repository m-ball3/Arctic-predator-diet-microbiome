# ------------------------------------------------------------------
# COMPARES 16S COMBINED ERROR RATES WITH BATCHED ERROR RATES
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Sets up the Environment and Loads in data
# ------------------------------------------------------------------
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

sessionInfo()
library(phyloseq);  packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2);   packageVersion("ggplot2")
library(tidyverse)
library(dplyr)
library(vegan)
library(Rfast)
library(ALDEx2)

# ------------------------------------------------------------------
# Load batched and combined phyloseq objects
# IMPORTANT: ASV labels may not match across objects, so we will align by sequence.
# ------------------------------------------------------------------

load("./Scripts/16s/rdata/16s.ps.RData")
ps.batched <- ps.filt          # batched DADA2 run

load("./Scripts/16s/rdata/16s.ps-test.RData")
ps.combined <- ps.filt         # combined DADA2 run

# Make sure both have refseq (ASV sequences) stored
# (If already done earlier, you can skip this; otherwise attach sequences here)
# Example if taxa_names are currently sequences:
# dna_batched  <- DNAStringSet(taxa_names(ps.batched))
# names(dna_batched) <- taxa_names(ps.batched)
# refseq(ps.batched) <- dna_batched
# taxa_names(ps.batched) <- paste0("ASV", seq_len(ntaxa(ps.batched)))
#
# dna_combined <- DNAStringSet(taxa_names(ps.combined))
# names(dna_combined) <- taxa_names(ps.combined)
# refseq(ps.combined) <- dna_combined
# taxa_names(ps.combined) <- paste0("ASV", seq_len(ntaxa(ps.combined)))

# ------------------------------------------------------------------
# Align taxa across methods by SEQUENCE, not ASV number
# ------------------------------------------------------------------

# Extract ASV maps: ASV ID, sequence, taxonomy
asv_map_batched <- data.frame(
  ASV      = taxa_names(ps.batched),
  sequence = as.character(refseq(ps.batched)),
  as.data.frame(tax_table(ps.batched)),
  stringsAsFactors = FALSE
)

asv_map_combined <- data.frame(
  ASV      = taxa_names(ps.combined),
  sequence = as.character(refseq(ps.combined)),
  as.data.frame(tax_table(ps.combined)),
  stringsAsFactors = FALSE
)

# Build a global sequence-based alignment: one row per unique sequence
asv_alignment <- full_join(
  asv_map_batched %>% dplyr::select(ASV_batched = ASV, sequence),
  asv_map_combined %>% dplyr::select(ASV_combined = ASV, sequence),
  by = "sequence"
)

# Add global ASV IDs (consistent across both methods)
unique_seqs <- unique(asv_alignment$sequence)
global_ASV  <- paste0("gASV", seq_along(unique_seqs))
seq_to_global <- data.frame(
  sequence  = unique_seqs,
  global_ASV = global_ASV,
  stringsAsFactors = FALSE
)

# Map global ASV IDs back onto both phyloseq objects
asv_alignment <- asv_alignment %>%
  left_join(seq_to_global, by = "sequence")

# Create lookup: sequence -> global_ASV
global_lookup <- setNames(seq_to_global$global_ASV, seq_to_global$sequence)

# Re-label taxa_names in ps.batched and ps.combined by global_ASV
taxa_names(ps.batched)  <- global_lookup[as.character(refseq(ps.batched))]
taxa_names(ps.combined) <- global_lookup[as.character(refseq(ps.combined))]

# Now taxa_names(ps.batched) and taxa_names(ps.combined) represent the same sequences with the same IDs.

# ------------------------------------------------------------------
# Restrict to overlapping samples and non-zero counts
# ------------------------------------------------------------------

batched.samps   <- sample_names(ps.batched)
combined.samps  <- sample_names(ps.combined)
sam.both        <- intersect(batched.samps, combined.samps)

ps.batched.abs  <- prune_samples(sam.both, ps.batched)
ps.combined.abs <- prune_samples(sam.both, ps.combined)

keep_samples <- sam.both[
  sample_sums(ps.batched.abs)[sam.both]  > 0 &
    sample_sums(ps.combined.abs)[sam.both] > 0
]

ps.batched.abs  <- prune_samples(keep_samples, ps.batched.abs)
ps.combined.abs <- prune_samples(keep_samples, ps.combined.abs)

# ------------------------------------------------------------------
# Relative abundance, combined df, and visualization
# ------------------------------------------------------------------

rel.batched  <- transform_sample_counts(ps.batched.abs,  function(x) x / sum(x))
rel.combined <- transform_sample_counts(ps.combined.abs, function(x) x / sum(x))

df_batched   <- psmelt(rel.batched)
df_combined  <- psmelt(rel.combined)
df_batched$type   <- "batched"
df_combined$type  <- "combined"

# Ensure Species column exists consistently
# (If tax_glom to Species was used earlier, Species should be present in tax_table)
combined_df <- bind_rows(df_batched, df_combined)

equal.one <- combined_df %>%
  group_by(Sample, type) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  print(n = 50)

sample_data(rel.batched)$type   <- "batched"
sample_data(rel.combined)$type  <- "combined"
sample_data(ps.batched.abs)$type   <- "batched"
sample_data(ps.combined.abs)$type  <- "combined"

summary(sample_sums(rel.batched))
summary(sample_sums(rel.combined))

# Plot comparison
rel.plot <- ggplot(combined_df, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_col(position = "stack") +
  facet_wrap(. ~ type, ncol = 1, strip.position = "right") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")

rel.plot

#------------------------------------------------------------------------------------
# ALPHA DIVERSITY
#------------------------------------------------------------------------------------

alpha_batched  <- estimate_richness(ps.batched.abs,  measures = c("Shannon", "Observed"))
alpha_combined <- estimate_richness(ps.combined.abs, measures = c("Shannon", "Observed"))

alpha_batched$Sample  <- rownames(alpha_batched)
alpha_combined$Sample <- rownames(alpha_combined)
alpha_batched$method  <- "batched"
alpha_combined$method <- "combined"

alpha_all <- bind_rows(alpha_batched, alpha_combined)

ggplot(alpha_all, aes(x = method, y = Shannon, group = Sample)) +
  geom_point() +
  geom_line(alpha = 0.3) +
  theme_minimal()

alpha_all$method <- factor(alpha_all$method, levels = c("batched", "combined"))

# Global test
kw_shannon <- kruskal.test(Shannon ~ method, data = alpha_all)
kw_shannon

# Paired test
alpha_wide <- alpha_all %>%
  dplyr::select(Sample, method, Shannon) %>%
  pivot_wider(names_from = method, values_from = Shannon)

wilcox_shannon <- wilcox.test(alpha_wide$batched, alpha_wide$combined, paired = TRUE)
wilcox_shannon

alpha_wide %>%
  mutate(diff = batched - combined) %>%
  summarise(median_diff = median(diff), mean_diff = mean(diff))

#------------------------------------------------------------------------------------
# BETA DIVERSITY (robust Aitchison) with paired design
#------------------------------------------------------------------------------------

# 1. Extract matrices: samples x global_ASV
otu_batched <- as(otu_table(rel.batched), "matrix")
otu_combined <- as(otu_table(rel.combined), "matrix")

if (taxa_are_rows(rel.batched))   otu_batched  <- t(otu_batched)
if (taxa_are_rows(rel.combined))  otu_combined <- t(otu_combined)

# 2. Align taxa (global_ASV) across methods by columns
all_taxa <- union(colnames(otu_batched), colnames(otu_combined))

otu_batched_aligned <- matrix(
  0,
  nrow = nrow(otu_batched),
  ncol = length(all_taxa),
  dimnames = list(rownames(otu_batched), all_taxa)
)

otu_comb_aligned <- matrix(
  0,
  nrow = nrow(otu_combined),
  ncol = length(all_taxa),
  dimnames = list(rownames(otu_combined), all_taxa)
)

idx_batched  <- match(colnames(otu_batched),  all_taxa)
idx_combined <- match(colnames(otu_combined), all_taxa)

otu_batched_aligned[, idx_batched] <- otu_batched
otu_comb_aligned[, idx_combined]   <- otu_combined

# 3. Stack rows: each sample-method combination is a row
otu_stack <- rbind(otu_batched_aligned, otu_comb_aligned)

samples_batched  <- rownames(otu_batched_aligned)
samples_combined <- rownames(otu_comb_aligned)

meta_stack <- data.frame(
  SampleID = c(samples_batched, samples_combined),
  method   = c(rep("batched",  length(samples_batched)),
               rep("combined", length(samples_combined)))
)

stopifnot(all(rownames(otu_stack) == meta_stack$SampleID))

# 4. Robust Aitchison distance
robust_aitch_dist <- vegdist(otu_stack, method = "robust.aitchison")

# Global effect of method
permanova_aitch_global <- adonis2(
  robust_aitch_dist ~ method,
  data        = meta_stack,
  permutations = 999
)

# Paired effect (block by SampleID)
permanova_aitch_paired <- adonis2(
  robust_aitch_dist ~ method,
  data        = meta_stack,
  strata      = meta_stack$SampleID,
  permutations = 999
)

print(permanova_aitch_global)
print(permanova_aitch_paired)

#------------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE (ALDEx2) aligned by global_ASV
#------------------------------------------------------------------------------------

counts_batched  <- as(otu_table(ps.batched.abs), "matrix")
counts_combined <- as(otu_table(ps.combined.abs), "matrix")

if (taxa_are_rows(ps.batched.abs))   counts_batched  <- t(counts_batched)
if (taxa_are_rows(ps.combined.abs))  counts_combined <- t(counts_combined)

# Align taxa (rows = global_ASV) across methods
all_taxa_counts <- union(rownames(counts_batched), rownames(counts_combined))

counts_batched_aligned <- matrix(
  0,
  nrow = length(all_taxa_counts),
  ncol = ncol(counts_batched),
  dimnames = list(all_taxa_counts, colnames(counts_batched))
)

counts_combined_aligned <- matrix(
  0,
  nrow = length(all_taxa_counts),
  ncol = ncol(counts_combined),
  dimnames = list(all_taxa_counts, colnames(counts_combined))
)

idx_batched_counts  <- match(rownames(counts_batched),  all_taxa_counts)
idx_combined_counts <- match(rownames(counts_combined), all_taxa_counts)

counts_batched_aligned[idx_batched_counts, ]  <- counts_batched
counts_combined_aligned[idx_combined_counts, ] <- counts_combined

# Stack columns: taxa x (sample × method)
counts_stack <- cbind(counts_batched_aligned, counts_combined_aligned)

group <- c(rep("batched",  ncol(counts_batched_aligned)),
           rep("combined", ncol(counts_combined_aligned)))

sample_ids <- c(colnames(counts_batched_aligned),
                colnames(counts_combined_aligned))

sample_ids <- paste(sample_ids, group, sep = "_")
colnames(counts_stack) <- sample_ids

aldex_out <- aldex(
  counts_stack,
  group,
  test   = "t",
  effect = TRUE,
  denom  = "all"   # avoid IQLR issues in sparse data
)

head(aldex_out)

sig_taxa <- subset(aldex_out, we.eBH < 0.05 & abs(effect) > 1)
head(sig_taxa)

# Annotate with taxonomy (global_ASV)
tax_batched_df <- as.data.frame(tax_table(ps.batched.abs))
tax_batched_df$ASV <- rownames(tax_batched_df)

sig_taxa_df <- as.data.frame(sig_taxa)
sig_taxa_df$ASV <- rownames(sig_taxa_df)

sig_annot_batched <- merge(sig_taxa_df, tax_batched_df, by = "ASV")
head(sig_annot_batched)

tax_combined_df <- as.data.frame(tax_table(ps.combined.abs))
tax_combined_df$ASV <- rownames(tax_combined_df)

sig_annot_combined <- merge(sig_taxa_df, tax_combined_df, by = "ASV")
head(sig_annot_combined)

#------------------------------------------------------------------------------------
# Per-sample composition differences (using global_ASV)
#------------------------------------------------------------------------------------

otu_batched_rel  <- as(otu_table(rel.batched), "matrix")
otu_combined_rel <- as(otu_table(rel.combined), "matrix")

if (taxa_are_rows(rel.batched))   otu_batched_rel  <- t(otu_batched_rel)
if (taxa_are_rows(rel.combined))  otu_combined_rel <- t(otu_combined_rel)

# Align taxa (rows = global_ASV)
all_taxa_rel <- union(rownames(otu_batched_rel), rownames(otu_combined_rel))

batched_rel_aligned <- matrix(
  0,
  nrow = length(all_taxa_rel),
  ncol = ncol(otu_batched_rel),
  dimnames = list(all_taxa_rel, colnames(otu_batched_rel))
)

combined_rel_aligned <- matrix(
  0,
  nrow = length(all_taxa_rel),
  ncol = ncol(otu_combined_rel),
  dimnames = list(all_taxa_rel, colnames(otu_combined_rel))
)

idx_batched_rel  <- match(rownames(otu_batched_rel),  all_taxa_rel)
idx_combined_rel <- match(rownames(otu_combined_rel), all_taxa_rel)

batched_rel_aligned[idx_batched_rel, ]   <- otu_batched_rel
combined_rel_aligned[idx_combined_rel, ] <- otu_combined_rel

common_samples <- intersect(colnames(batched_rel_aligned), colnames(combined_rel_aligned))

equal_vec <- function(x, y, tol = 1e-8) {
  all(abs(x - y) < tol)
}

sample_equal <- sapply(common_samples, function(s) {
  x <- batched_rel_aligned[, s]
  y <- combined_rel_aligned[, s]
  equal_vec(x, y)
})

samples_same <- names(sample_equal)[sample_equal]
samples_diff <- names(sample_equal)[!sample_equal]

samples_same
samples_diff

max_diff <- sapply(common_samples, function(s) {
  x <- batched_rel_aligned[, s]
  y <- combined_rel_aligned[, s]
  max(abs(x - y))
})

diff_table <- data.frame(
  Sample   = common_samples,
  equal    = sample_equal,
  max_diff = max_diff
)

diff_table

ggplot(diff_table, aes(x = Sample, y = max_diff)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

# Build ASV-level difference table for samples with max_diff >= 0.25
samples_to_keep <- diff_table %>%
  filter(max_diff >= 0.25) %>%
  pull(Sample)

asv_diff_all <- lapply(samples_to_keep, function(s) {
  x <- batched_rel_aligned[, s]
  y <- combined_rel_aligned[, s]
  
  diff_vec <- x - y
  
  asv_diff <- data.frame(
    Taxon        = rownames(batched_rel_aligned),
    rel_batched  = x,
    rel_combined = y,
    diff         = diff_vec,
    stringsAsFactors = FALSE
  )
  
  asv_diff_nonzero <- asv_diff %>%
    filter(abs(diff) > 1e-8)
  
  asv_diff_nonzero$Sample <- s
  
  asv_diff_nonzero
})

asv_diff_table <- bind_rows(asv_diff_all)

asv_diff_table <- asv_diff_table %>%
  dplyr::rename(SampleID = "Taxon")
asv_diff_table <- asv_diff_table %>%
  dplyr::rename(Taxon = "Sample")

# Taxonomy (Genus/Species) by global_ASV
tax_batched_simple <- as.data.frame(tax_table(ps.batched.abs))
tax_batched_simple$Taxon <- rownames(tax_batched_simple)
tax_batched_simple <- tax_batched_simple %>%
  dplyr::select(Taxon, Genus, Species)

tax_combined_simple <- as.data.frame(tax_table(ps.combined.abs))
tax_combined_simple$Taxon <- rownames(tax_combined_simple)
tax_combined_simple <- tax_combined_simple %>%
  dplyr::select(Taxon, Genus, Species)

asv_diff_tax <- asv_diff_table %>%
  left_join(tax_batched_simple,  by = "Taxon", suffix = c("", "_batched")) %>%
  left_join(tax_combined_simple, by = "Taxon", suffix = c("", "_combined")) %>%
  dplyr::rename(Genus_batched = "Genus")%>%
  dplyr::rename(Species_batched = "Species")%>%
  dplyr::mutate(
    same_genus   = Genus_batched   == Genus_combined,
    same_species = Species_batched == Species_combined
  )


# Write out
getwd()
write.csv(asv_diff_tax, "./Deliverables/16S/asv_diff_tax-learnerrorrates.csv", row.names = FALSE)
