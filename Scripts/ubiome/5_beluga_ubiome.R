# ------------------------------------------------------------------
# CREATES PHYLOSEQ OBJS FOR BELUGA SAMPLES - INTENDED FOR FREDDY'S SAMPLE COMPARISONS
# THIS IS AN OPTIONAL FIFTH STEP AFTER DADA2
## 1_rownames-match_ubiome.R and 
## 2_replicates_contaminated_ubiome.R 3_ampbias_ubiome.R
## and 4_phyloseq_ubiome.R must be run before this
# ------------------------------------------------------------------

# Sets up environment
library(phyloseq)

# loads in phyloseq objects
load("ubiome.ps-major.Rdata")
load("ubiome.ps.Rdata")


## Adds a new column to the ps tax tables
# # Extract the taxonomy table as a matrix
# tax1 <- tax_table(ps.raw)
# tax2 <- tax_table(ps.major)
# 
# # Create the new column: "Phylum: Family"
# phylum_family1 <- paste(tax1[, "Phylum"], tax1[, "Family"], sep = ": ")
# phylum_family2 <- paste(tax2[, "Phylum"], tax2[, "Family"], sep = ": ")
# 
# # 3. Add it back to the taxonomy table
# tax_table(ps.raw) <- cbind(tax1, `Phylum: Family` = phylum_family1)
# tax_table(ps.major) <- cbind(tax2, `Phylum: Family` = phylum_family2)

# Filters samdf for only beluga samples
samdf.beluga <- samdf %>%
  filter(Predator == "beluga whale")

# RAW ---------------------------------------------------------------------------------------

# Filters for only beluga samples
ps.beluga <- subset_samples(ps.raw, Predator == "beluga whale")

## MERGE TO SPECIES HERE (TAX GLOM)
ps.beluga = tax_glom(ps.beluga, "Family", NArm = FALSE)

# transforms to proportional abundance
ps.beluga.prop <- transform_sample_counts(ps.beluga, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Plots with WADE IDs - MAJOR
beluga.prop.plot <- plot_bar(ps.beluga.prop, fill="Phylum")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
beluga.prop.plot

# MAJOR--------------------------------------------------------------------------------------
# Filters for only beluga samples
ps.beluga.major <- subset_samples(ps.major, Predator == "beluga whale")

## MERGE TO SPECIES HERE (TAX GLOM)
ps.belugamajor = tax_glom(ps.beluga.major, "Family", NArm = FALSE)


# transforms to proportional abundance
ps.beluga.major.prop <- transform_sample_counts(ps.beluga.major, function(x) {
  x_rel <- x / sum(x)
  x_rel[is.nan(x_rel)] <- 0
  return(x_rel)
})

# Plots with WADE IDs - MAJOR
beluga.prop.plot.major <- plot_bar(ps.beluga.major.prop, fill="Phylum")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
beluga.prop.plot.major

#saves plots 
ggsave("Deliverables/ubiome/ubiome-beluga-majorphylum.png", plot = beluga.prop.plot.major, width = 30, height = 8, units = "in", dpi = 300)
ggsave("Deliverables/ubiome/ubiome-beluga-allphylum.png", plot = beluga.prop.plot, width = 30, height = 8, units = "in", dpi = 300)


# ------------------------------------------------------------------
# TABLES
# ------------------------------------------------------------------

# CREATES ABSOLUTE SAMPLES X SPECIES TABLE 
otu.abs <- as.data.frame(otu_table(ps.beluga.major))
colnames(otu.abs) <- as.data.frame(tax_table(ps.beluga.major))$Family

# Checks for differences = WADE-003-178
setdiff(rownames(samdf.beluga), rownames(otu.abs))
setdiff(rownames(otu.abs), rownames(samdf.beluga))

## Adds ADFG Sample ID as a column
otu.abs$Specimen.ID <- samdf.beluga$Specimen.ID # error here

## Moves ADFG_SampleID to the first column
otu.abs <- otu.abs[, c(ncol(otu.abs), 1:(ncol(otu.abs)-1))]

# CREATES RELATIVE SAMPLES X SPECIES TABLE
otu.prop <- as.data.frame(otu_table(ps.beluga.major.prop))
colnames(otu.prop) <- as.data.frame(tax_table(ps.beluga.major.prop))$Family

## Adds ADFG Sample ID as a column (do NOT set as row names if not unique)
otu.prop$Specimen.ID <- samdf.beluga$Specimen.ID

## Moves ADFG_SampleID to the first column
otu.prop <- otu.prop[, c(ncol(otu.prop), 1:(ncol(otu.prop)-1))]

# Changes NaN to 0
#otu.prop[is.na(otu.prop)] <- 0

# Rounds to three decimal places
is.num <- sapply(otu.prop, is.numeric)
otu.prop[is.num] <- lapply(otu.prop[is.num], round, 3)

# Writes to CSV
write.csv(otu.abs, "./Deliverables/ubiome/ubiome_absolute_speciesxsamples-BELUGAMAJOR.csv", row.names = TRUE)
write.csv(otu.prop, "./Deliverables/ubiome/ubiome_relative_speciesxsamples-BELUGAMAJOR.csv", row.names = TRUE)

