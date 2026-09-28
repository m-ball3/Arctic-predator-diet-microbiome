# ------------------------------------------------------------------------------
# MERGES 12S REGIONAL PHYLOSEQ OBJECTS
# 1_rownames_match_ubiome.R, 2_decontam_ubiome.R, 
# 3_replicates_ubiome.R, 4_phyloseq_ubiome.R should all be run before this
# ------------------------------------------------------------------------------

# loads data
# Loads dada2 output
# load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-COOKINLET.Rdata")
# load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-SBERING.Rdata")
# load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-ARCTIC.Rdata")

load("./Scripts/12s/rdata/ps.12s.regions.raw.Rdata")

metadata.arctic <- data.frame(sample_data(ps.12s.arctic))

# sanity checks - we want sequences here
head(phyloseq::taxa_names(ps.12s.cook))
head(phyloseq::taxa_names(ps.12s.sbering))
head(phyloseq::taxa_names(ps.12s.arctic))

# sanity checks - sample names
head(phyloseq::sample_names(ps.12s.cook))
head(phyloseq::sample_names(ps.12s.sbering))
head(phyloseq::sample_names(ps.12s.arctic))

identical(
  sample_names(ps.12s.cook),
  rownames(as.data.frame(sample_data(ps.12s.cook)))
)

identical(
  sample_names(ps.12s.sbering),
  rownames(as.data.frame(sample_data(ps.12s.sbering)))
)

identical(
  sample_names(ps.12s.arctic),
  rownames(as.data.frame(sample_data(ps.12s.arctic)))
)

# merges phyloseq objects
ps.12s.merged.asv <- merge_phyloseq(
  ps.12s.cook,
  ps.12s.sbering,
  ps.12s.arctic
)


# sanity checks for correct merging----------------------------------------------------------------

# Are the expected number of samples present after merging?
n.expected <- phyloseq::nsamples(ps.12s.cook) +
  phyloseq::nsamples(ps.12s.sbering) +
  phyloseq::nsamples(ps.12s.arctic)

n.merged <- phyloseq::nsamples(ps.12s.merged.asv)

c(
  cook = phyloseq::nsamples(ps.12s.cook),
  sbering = phyloseq::nsamples(ps.12s.sbering),
  arctic = phyloseq::nsamples(ps.12s.arctic),
  expected_total = n.expected,
  merged_total = n.merged
)

# Are all regional samples present in the merged object?
all.original.samples <- c(
  phyloseq::sample_names(ps.12s.cook),
  phyloseq::sample_names(ps.12s.sbering),
  phyloseq::sample_names(ps.12s.arctic)
)


all(all.original.samples %in% phyloseq::sample_names(ps.12s.merged.asv))

# Which original samples, if any, are missing?
setdiff(
  all.original.samples,
  phyloseq::sample_names(ps.12s.merged.asv)
)

# Are there unexpected sample IDs in the merged object?
setdiff(
  phyloseq::sample_names(ps.12s.merged.asv),
  all.original.samples
)

# did the ASVs merge correctly?
taxa.expected <- union(
  union(
    taxa_names(ps.12s.cook),
    taxa_names(ps.12s.sbering)
  ),
  taxa_names(ps.12s.arctic)
)

c(
  expected_unique_ASVs = length(taxa.expected),
  merged_ASVs = phyloseq::ntaxa(ps.12s.merged.asv)
)

setdiff(taxa.expected, taxa_names(ps.12s.merged.asv))
setdiff(taxa_names(ps.12s.merged.asv), taxa.expected)

# are counts within samples the same after merging?
check_sample_counts <- function(source_ps, merged_ps) {
  
  source.otu <- as(phyloseq::otu_table(source_ps), "matrix")
  merged.otu <- as(phyloseq::otu_table(merged_ps), "matrix")
  
  # Ensure matrices are samples × ASVs
  if (phyloseq::taxa_are_rows(source_ps)) {
    source.otu <- t(source.otu)
  }
  
  if (phyloseq::taxa_are_rows(merged_ps)) {
    merged.otu <- t(merged.otu)
  }
  
  # Subset merged counts to the exact samples and ASVs in source object
  merged.subset <- merged.otu[
    rownames(source.otu),
    colnames(source.otu),
    drop = FALSE
  ]
  
  # Confirm matching dimensions and names before comparing values
  stopifnot(
    identical(rownames(source.otu), rownames(merged.subset)),
    identical(colnames(source.otu), colnames(merged.subset))
  )
  
  # Compare count values; NA-safe
  all(
    source.otu == merged.subset,
    na.rm = FALSE
  )
}
check_sample_counts(ps.12s.cook, ps.12s.merged.asv)

check_sample_counts(ps.12s.sbering, ps.12s.merged.asv)

check_sample_counts(ps.12s.arctic, ps.12s.merged.asv)

# IF ABOVE IS FALSE, -----------------------------------------------------------------------
compare_sample_counts <- function(source_ps, merged_ps) {
  
  source.otu <- as(otu_table(source_ps), "matrix")
  merged.otu <- as(otu_table(merged_ps), "matrix")
  
  if (taxa_are_rows(source_ps)) source.otu <- t(source.otu)
  if (taxa_are_rows(merged_ps)) merged.otu <- t(merged.otu)
  
  merged.subset <- merged.otu[
    rownames(source.otu),
    colnames(source.otu),
    drop = FALSE
  ]
  
  mismatch <- source.otu != merged.subset
  
  list(
    n_mismatched_cells = sum(mismatch),
    mismatches = which(mismatch, arr.ind = TRUE)
  )
}

# we want 0 mismatches
compare_sample_counts(ps.12s.cook, ps.12s.merged.asv)

compare_sample_counts(ps.12s.sbering, ps.12s.merged.asv)

compare_sample_counts(ps.12s.arctic, ps.12s.merged.asv)


#-------------------------------------------------------------------------------------


# Filters out anything not in Actinopteri
ps.12s.merged.asv <- phyloseq::subset_taxa(ps.12s.merged.asv , Class == "Actinopteri")
phyloseq::nsamples(ps.12s.merged.asv)

# Remove samples with total abundance < 100 <- DO I MEAN READS HERE???
ps.12s.merged <- phyloseq::prune_samples(phyloseq::sample_sums(ps.12s.merged.asv) >= 100, ps.12s.merged.asv)
phyloseq::sample_sums(ps.12s.merged)
phyloseq::nsamples(ps.12s.merged)

# Then, if appropriate:
ps.12s.merged <- tax_glom(
  ps.12s.merged,
  taxrank = "Species",
  NArm = FALSE
)

metadata.12s.merged <- data.frame(sample_data(ps.12s.merged))


# Saves phyloseq obj per region (RAW)
save(ps.12s.merged, metadata.12s.merged, file = "./Scripts/12s/rdata/ps.12s.merged.Rdata")
