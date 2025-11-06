## Table append

# Assuming your three tables loaded from RData file:
# genus.species (column: genus, species)
# taxa (multiple taxonomic columns e.g. Kingdom, Phylum, etc.)
# species (taxa table with appended species1 column)

# Convert matrices to data frames so columns can be added
taxa_df <- as.data.frame(taxa)
genus_species_df <- as.data.frame(genus.species)
species_df <- as.data.frame(species)

# Create new column with combined genus and species binomial (space-separated)
genus_species_df$species_binomial <- ifelse(
  is.na(genus_species_df$Genus) | is.na(genus_species_df$Species),
  NA,
  paste(genus_species_df$Genus, genus_species_df$Species, sep = " ")
)

# Optionally remove columns if not needed
genus_species_df$Genus <- NULL
genus_species_df$Species <- NULL
species_df$Kingdom <- NULL
species_df$Phylum <- NULL
species_df$Class <- NULL
species_df$Order <- NULL
species_df$Family <- NULL
species_df$Genus <- NULL
species_df$Species <- NULL

# Add Sequence as a cNULL# Add Sequence as a column with rownames
taxa_df$Sequence <- rownames(taxa_df)
genus_species_df$Sequence <- rownames(genus_species_df)
species_df$Sequence <- rownames(species_df)

# Rename columns to indicate method/source
colnames(taxa_df)[colnames(taxa_df) != "Sequence"] <- paste0(colnames(taxa_df)[colnames(taxa_df) != "Sequence"], "-assigntaxonomy")
colnames(genus_species_df)[colnames(genus_species_df) != "Sequence"] <- paste0(colnames(genus_species_df)[colnames(genus_species_df) != "Sequence"], "-assignspecies")
colnames(species_df)[colnames(species_df) != "Sequence"] <- paste0(colnames(species_df)[colnames(species_df) != "Sequence"], "-addspecies")

# Merge by Sequence column
merged_table <- merge(taxa_df, genus_species_df, by = "Sequence", all = TRUE)
merged_table <- merge(merged_table, species_df, by = "Sequence", all = TRUE)

# Set Sequence back as rownames then remove Sequence column
rownames(merged_table) <- merged_table$Sequence
merged_table$Sequence <- NULL

class(merged_table)

save(seqtab.nochim, freq.nochim, track, taxa, genus.species, species, merged_table, file = "WADE003-arcticpred_dada2_QAQC_12SP1_output-addSpecies-130;30-3.Rdata")
