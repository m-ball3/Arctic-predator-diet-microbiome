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
load("./Scripts/ubiome/rdata/rownames-match_ubiome.RData")

# ------------------------------------------------------------------------------
# Contamination
# ------------------------------------------------------------------------------
## NOW THINK ABOUT NEGATIVE CONTAMINATION!!!!!

# ps for negs
samdfneg <- samdf %>% 
  filter(grepl("neg", rownames(.)))

seqtab.negtest <- as.data.frame(seqtab.nochim) %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample, names_to = "seq", values_to = "nReads") %>% 
  mutate(sample_type = case_when(grepl("mock", sample)~"mock",
                                 grepl("neg", sample)~"neg",
                                 TRUE~"sample")) %>% 
  group_by(seq, sample_type) %>% 
  summarise(nReadsTotal = sum(nReads)) %>% 
  ungroup() %>% group_by(seq) %>% 
  mutate(seqNumber = cur_group_id()) %>% 
  filter(!any(sample_type == "neg" & nReadsTotal == 0)) %>% 
  arrange(desc(nReadsTotal)) %>% 
  filter(sample_type != "mock") %>% 
  filter(!any(sample_type == "sample" & nReadsTotal == 0))

# 
# neg_contaminated <- 
#   #some vector here to use to filter for the contaminated samples from the shipment 5 extractions
# 
# # need to change this syntax because I don't want to remove them, 
#   ## instead, I want to delete the reads that appear in both
#   ### MAYBE THIS SHOULD HAPPEN BEFORE THE REPLICATE STEPS???
# samdf_filt <- samdf_filt[!rownames(samdf_filt) %in% neg_contaminated, ]
# 
# seqtab.nochim_filt <- seqtab.nochim_filt[!rownames(seqtab.nochim_filt) %in% neg_contaminated, ]
# 

# ------------------------------------------------------------------------------
# Replicates
# ------------------------------------------------------------------------------

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
plot_bar(ps.replicates, x = "LabID", fill = "Genus") +
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
  "WADE-003-151",
  "WADE-003-128",
  "WADE-003-170",
  "WADE-003-111",
  "WADE-003-146",
  "WADE-003-123",
  "WADE-003-216",
  "WADE-003-201",
  "WADE-003-215"
)

samdf <- samdf[!rownames(samdf) %in% replicate_to_remove, ]
seqtab.nochim <- seqtab.nochim[!rownames(seqtab.nochim) %in% replicate_to_remove, ]

# Recreates phyloseq object without unwanted replicates
ps.norep <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
  sample_data(samdf), 
  tax_table(taxam)
) %>%
  subset_samples(!LabID %in% replicate_to_remove)

sample_names(ps.norep)
nsamples(ps.norep)


# REMOVES TISSUE SAMPLE (NOT AN SRKW FECAL SAMPLE)
# DEALS WITH REPLICATES
## Visually compare their absolute and proportional abundance in a separate bar plot. 
## You can also run NMDS on them to see if they cluster. 
## This will tell us a little bit about sequencing/subsampling bias on our data.

getwd()
save(samdf, seqtab.nochim, taxam, track, out, freq.nochim, ps.raw, ps.norep, file = "./Scripts/ubiome/rdata/replicates-contaminated_ubiome.RData")
