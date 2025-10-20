# Assuming your data is in a dataframe called df
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load the object, don't assign
load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-130trunc-ADFGnotes.Rdata")
# Removes file extensions from OTU table names
rownames(track) <- gsub("-MFU_S\\d+", "", rownames(track))

# List all objects in the environment
ls()

# Transformes track to df
track.df <- as.data.frame(track)

# Move sample names from rownames to a column
track.12s <- tibble::rownames_to_column(track.df, var = "Sample")

# Pivot longer to get a tidy dataframe for ggplot2
track_long <- pivot_longer(
  track.12s,
  cols = -Sample,
  names_to = "Step",
  values_to = "Reads"
)

# Ensure your data is in long format as shown before
# Use your wide-format table to identify qualifying samples
high_samples <- track.12s$Sample[track.12s$nonchim >= 50000]
low_samples  <- track.12s$Sample[track.12s$nonchim < 50000]

# Filter long-format data by sample sets
track_high <- track_long[track_long$Sample %in% high_samples, ]
track_low  <- track_long[track_long$Sample %in% low_samples, ]

# Plot for samples >= 35000 final reads
p1 <- ggplot(track_high, aes(x = Sample, y = Reads, fill = Step)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Number of Reads", fill = "Step", title = "Samples with 50000+ Nonchim Reads")

# Plot for samples < 35000 final reads
p2 <- ggplot(track_low, aes(x = Sample, y = Reads, fill = Step)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Number of Reads", fill = "Step", title = "Samples with <50000 Nonchim Reads")

# Plot for nochim read count
# "High" and "Low" subsets for nonchim plot data frames
track_high_nchim <- track.12s[track.12s$Sample %in% high_samples, ]
track_low_nchim  <- track.12s[track.12s$Sample %in% low_samples, ]

# Plot for high nonchim
p_high_nchim <- ggplot(track_high_nchim, aes(x = Sample, y = nonchim)) +
  geom_bar(stat = "identity", fill = "hotpink") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Nonchim Reads", title = "Samples with 50000+ Nonchim Reads (Nonchim Only)")

# Plot for low nonchim
p_low_nchim <- ggplot(track_low_nchim, aes(x = Sample, y = nonchim)) +
  geom_bar(stat = "identity", fill = "hotpink") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample", y = "Nonchim Reads", title = "Samples with <50000 Nonchim Reads (Nonchim Only)")

# Print your plots (using patchwork)
p1/p2 
p_high_nchim / p_low_nchim



