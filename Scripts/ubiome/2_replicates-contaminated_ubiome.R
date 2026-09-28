# ------------------------------------------------------------------
# DEALS WITH TECHNICAL REPLICATES AND ANY CONTAMINATED SAMPLES
# THIS IS THE SECOND STEP AFTER DADA2
## 1_rownames-match_ubiome.R must be run before this
## if there are contaminated samples,
## 1.5_decontam_ubiome.R must be run before this
# ------------------------------------------------------------------

## Sets up the Environment and Loads Libraries
library(tidyverse)
library(dplyr)
library(tibble)
library(stringr)
library(lubridate)
library(Biostrings)
library(openxlsx)

# Loads in output from DADA2 + filtered seqtab.nochim and samdf from rownames-match.r
load("./Scripts/ubiome/rdata/rownames-match_ubiome.RData")


# # samdf for negs and mock
# samdfnegmock <- samdf %>%
#   filter(
#     grepl("neg",  rownames(.), ignore.case = TRUE) |
#       grepl("mock", rownames(.), ignore.case = TRUE)
#   )
# 
# # get sample names that contain 'neg' or 'mock'
# neg_mock_ids <- sample_names(ps.raw)[
#   grepl("neg",  sample_names(ps.raw), ignore.case = TRUE) |
#     grepl("mock", sample_names(ps.raw), ignore.case = TRUE)
# ]
# 
# # subset to those samples
# ps.negmock <- prune_samples(neg_mock_ids, ps.raw)
# 
# ### shorten ASV seq names, store sequences as reference
# dna <- Biostrings::DNAStringSet(taxa_names(ps.negmock))
# names(dna) <- taxa_names(ps.negmock)
# ps.negmock <- merge_phyloseq(ps.negmock, dna)
# taxa_names(ps.negmock) <- paste0("ASV", seq(ntaxa(ps.negmock)))
# 
# nsamples(ps.negmock)
# 
# # Filters out any Mammalia
# ps.negmock <- subset_taxa(ps.negmock, Class!="Mammalia")
# 
# # Remove samples with total abundance == 0
# ps.negmock <- prune_samples(sample_sums(ps.negmock) > 0, ps.negmock)
# 
# ## MERGE TO SPECIES HERE (TAX GLOM)
# ps.negmock = tax_glom(ps.negmock, "Class", NArm = FALSE)
# 
# ## Transforms NaN (0/0) to 0
# ps.negmock.prop <- transform_sample_counts(ps.negmock, function(x) {
#   x_rel <- x / sum(x)
#   x_rel[is.nan(x_rel)] <- 0
#   return(x_rel)
# })
# 
# # PLOTS ------------------------------------------------------------------------
# # Plots with WADE IDs
# rel.plot <- plot_bar(ps.negmock.prop, fill="Class")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
#         legend.position = 'right')
# rel.plot
# 
# # ggsave("Deliverables/ubiome/ubiome-negmock-Class.png", plot = rel.plot, width = 30, height = 30, units = "in", dpi = 300)
# 
# # TABLES -----------------------------------------------------------------------
# # CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
# otu.abs <- as.data.frame(otu_table(ps.negmock))
# colnames(otu.abs) <- as.data.frame(tax_table(ps.negmock))$Family
# 
# # Checks for differences = WADE-003-178
# setdiff(rownames(samdfnegmock), rownames(otu.abs))
# setdiff(rownames(otu.abs), rownames(samdfnegmock))
# 
# ## Adds ADFG Sample ID as a column
# otu.abs$Specimen.ID <- samdfnegmock$Specimen.ID # error here
# 
# ## Moves ADFG_SampleID to the first column
# otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]
# 
# # CREATES RELATIVE SAMPLES X SPECIES TABLE
# otu.prop <- as.data.frame(otu_table(ps.negmock.prop))
# colnames(otu.prop) <- as.data.frame(tax_table(ps.negmock.prop))$Family
# 
# ## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
# otu.prop$Specimen.ID <- samdfnegmock$Specimen.ID
# 
# ## Moves ADFG_SampleID to the first column
# otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]
# 
# # Changes NaN to 0
# #otu.prop[is.na(otu.prop)] <- 0
# 
# # Rounds to three decimal places
# is.num <- sapply(otu.prop, is.numeric)
# otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)
# 
# # Writes to CSV
# write.csv(otu.abs, "./Deliverables/ubiome/ubiome_absolute_speciesxsamples-negmock-FAM.csv", row.names = TRUE)
# write.csv(otu.prop, "./Deliverables/ubiome/ubiome_relative_speciesxsamples-negmock-FAM.csv", row.names = TRUE)
# 
# 
# 
# seqtab.negtest <- as.data.frame(seqtab.nochim) %>% 
#   rownames_to_column("sample") %>% 
#   pivot_longer(-sample, names_to = "seq", values_to = "nReads") %>% 
#   mutate(sample_type = case_when(grepl("mock", sample)~"mock",
#                                  grepl("neg", sample)~"neg",
#                                  TRUE~"sample")) %>% 
#   group_by(seq, sample_type) %>% 
#   summarise(nReadsTotal = sum(nReads)) %>% 
#   ungroup() %>% group_by(seq) %>% 
#   mutate(seqNumber = cur_group_id()) %>% 
#   filter(!any(sample_type == "neg" & nReadsTotal == 0)) %>% 
#   arrange(desc(nReadsTotal)) %>% 
#   filter(sample_type != "mock") %>% 
#   filter(!any(sample_type == "sample" & nReadsTotal == 0))
# 
# # 
# # neg_contaminated <- 
# #   #some vector here to use to filter for the contaminated samples from the shipment 5 extractions
# # 
# # # need to change this syntax because I don't want to remove them, 
# #   ## instead, I want to delete the reads that appear in both
# #   ### MAYBE THIS SHOULD HAPPEN BEFORE THE REPLICATE STEPS???
# # samdf_filt <- samdf_filt[!rownames(samdf_filt) %in% neg_contaminated, ]
# # 
# # seqtab.nochim_filt <- seqtab.nochim_filt[!rownames(seqtab.nochim_filt) %in% neg_contaminated, ]
# # 
# 
# # saves for output into decontamination step
# save(samdf, seqtab.nochim, taxam, track, out, freq.nochim, ps.raw, ps.negmock, ps.negmock.prop, samdfnegmock, seqtab.negtest, file = "./Scripts/ubiome/rdata/replicates-contaminated-negmock_ubiome.RData")


# ------------------------------------------------------------------------------
# Replicates
# ------------------------------------------------------------------------------

#loads in data after decontamination step
load("./Scripts/12s/rdata/ps.decontaminated.RData")
# Determines which samples are replicates by pulling the duplicated specimen IDs
replicates <- unique(samdf$Specimen.ID[duplicated(samdf$Specimen.ID)])

# Subset phyloseq object to keep only samples with duplicated Specimen.IDs
ps.replicates <- subset_samples(ps.decontaminated, Specimen.ID %in% replicates)

# Creates a df
replicates_df <- samdf %>%
  dplyr::filter(Specimen.ID %in% replicates) %>%
  dplyr::select(Specimen.ID, LabID)

# Creates a vector of read counts
reads <- phyloseq::sample_sums(ps.replicates)

# Turn read counts into a data frame with a LabID column
reads_df <- data.frame(
  LabID = names(reads),
  reads = as.numeric(reads),
  row.names = NULL
)

# Join read counts onto your replicate table
replicates_df <- replicates_df %>%
  dplyr::left_join(reads_df, by = "LabID")

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
phyloseq::nsamples(ps.norep)


# DEALS WITH REPLICATES
## Visually compare their absolute and proportional abundance in a separate bar plot. 
## You can also run NMDS on them to see if they cluster. 
## This will tell us a little bit about sequencing/subsampling bias on our data.

getwd()
save(samdf, seqtab.nochim, taxam, track, out, freq.nochim, ps.raw, ps.norep, file = "./Scripts/ubiome/rdata/replicates-contaminated_ubiome.RData")
