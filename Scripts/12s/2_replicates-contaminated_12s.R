# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES AND ANY CONTAMINATED SAMPLES
# THIS IS THE SECOND STEP AFTER DADA2
## 1_rownames-match_ubiome.R must be run before this
# ------------------------------------------------------------------

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(stringr)
library(lubridate)


# Loads in output from DADA2 + filtered seqtab.nochim and samdf from rownames-match.r
load("./Scripts/12s/rdata/rownames-match_12s.RData")


# Determines which samples are replicates by pulling the duplicated specimen IDs
replicates <- unique(samdf$Specimen.ID[duplicated(samdf$Specimen.ID)])

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
ps.replicates <- subset_samples(ps.raw, Specimen.ID %in% replicates)

# Creates a df
replicates_df <- samdf %>%
  filter(Specimen.ID %in% replicates) %>%
  select(Specimen.ID, LabID)

# Creates a vector of read counts
reads <- sample_sums(ps.replicates)

# Turn read counts into a data frame with a LabID column
reads_df <- data.frame(
  LabID = names(reads),
  reads = as.numeric(reads),
  row.names = NULL
)

# Join read counts onto your replicate table
replicates_df <- replicates_df %>%
  left_join(reads_df, by = "LabID")

# Create the stacked bar plot (absolute)
plot_bar(ps.replicates, x = "LabID", fill = "Species") +
  facet_wrap(~ Specimen.ID, ncol = 13, scales = "free_x", strip.position = "top") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Transforms read counts to relative abundance of each species 
## Transforms NaN (0/0) to 0
ps.replicates.rel <- transform_sample_counts(ps.replicates, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Create the stacked bar plot (relative)
plot_bar(ps.replicates.rel, x = "LabID", fill = "Genus") +
  facet_wrap(~ Specimen.ID, ncol = 13, scales = "free_x", strip.position = "top") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Ensures samples removed in filtering are removed from samdf
replicate_to_remove <- c(
  "WADE-003-122",
  "WADE-003-219",
  "WADE-003-180",
  "WADE-003-124",
  "WADE-003-220",
  "WADE-003-190"
)

samdf <- samdf[!rownames(samdf) %in% replicate_to_remove, ]
seqtab.nochim <- seqtab.nochim[!rownames(seqtab.nochim) %in% replicate_to_remove, ]

# Recreates phyloseq object without unwanted replicates
ps.norep <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
  sample_data(samdf), 
  tax_table(taxa)
) %>%
  subset_samples(!LabID %in% replicate_to_remove)

sample_names(ps.norep)
nsamples(ps.norep) # = 122


# all 12s negatives were clean- check with Amy about needed to run contamination here

## NOW THINK ABOUT NEGATIVE CONTAMINATION!!!!!

# NEED TO PULL JUST THE DATE THAT THE CONTAMINATED EXTRACTION OCCURRED
## THE DATE SHOULD BE IN LABDF
# neg_contaminated <- 
#   #some vector here to use to filter for the contaminated samples from the shipment 5 extractions
#   
#   # need to change this syntax because I don't want to remove them, 
#   ## instead, I want to delete the reads that appear in both
#   ### MAYBE THIS SHOULD HAPPEN BEFORE THE REPLICATE STEPS???
#   samdf_filt <- samdf_filt[!rownames(samdf_filt) %in% neg_contaminated, ]
# 
# seqtab.nochim_filt <- seqtab.nochim_filt[!rownames(seqtab.nochim_filt) %in% neg_contaminated, ]


getwd()
save(samdf, seqtab.nochim, taxa, track, out, freq.nochim, ps.raw, ps.norep, file = "./Scripts/12s/rdata/replicates-contaminated_12s.RData")
