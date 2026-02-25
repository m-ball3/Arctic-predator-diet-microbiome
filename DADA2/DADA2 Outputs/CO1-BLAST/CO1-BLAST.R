library(dada2)
library(tidyverse)
library(stringr)
library(ggplot2)
library(seqinr)
library(dplyr)
library(Biostrings)
library(ShortRead)
library(rentrez)
library(xml2)

# Loads dada2 output
load("DADA2/DADA2 Outputs/testDADA2_CO1_allseqs_012826.Rdata")

# Renames columns to standard; removes _#### 
colnames(compiled_taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
compiled_taxa <- compiled_taxa %>%
  as.data.frame() %>%
  mutate(across(
    everything(),
    ~ str_replace(., "_.*$", "")   # remove underscore and everything after
  )) %>%
  as.matrix()

### Follow up on unassigned sequences
start_blast <- Sys.time()
# Find sequences with NA below kingdom
unassigned_idx <- which(is.na(compiled_taxa[, "Genus"]))  # or any level you care about

# Extract sequence names
seqs_unassigned <- colnames(seqtab.nochim)[unassigned_idx]

# Convert to DNAStringSet for export
seqs_unassigned_fasta <- DNAStringSet(seqs_unassigned)
names(seqs_unassigned_fasta) <- seqs_unassigned  # keep sequences as IDs

writeXStringSet(seqs_unassigned_fasta, "DADA2/DADA2 Outputs/CO1-BLAST/Co1_unassigned_seqs.fasta")

###uploaded this file to BLAST, then downloaded hits as a csv

unassigned_hits <- read_csv("DADA2/DADA2 Outputs/CO1-BLAST/TVYG56JC014-Alignment-HitTable.csv", col_names = FALSE)

accessions <- unique(unassigned_hits$X2)

tax_summaries <- lapply(accessions, function(acc) {
  summary <- entrez_summary(db="nuccore", id=acc)
  list(
    accession = acc,
    title = summary$title,
    taxid = summary$taxid,
    organism = summary$organism
  )
})

tax_df <- do.call(rbind, lapply(tax_summaries, as.data.frame))

unassigned_taxa <- unassigned_hits %>% 
  left_join(tax_df, by = c("X2" = "accession")) %>% 
  group_by(X1) %>% 
  separate(organism, into = c("Genus", "Species", "extra"), remove = FALSE) %>% 
  mutate(nGenus = length(unique(Genus))) %>% 
  mutate(nSpecies = length(unique(Species))) %>% 
  mutate(Species = case_when(nSpecies>1~"sp.",
                             TRUE~Species)) %>% 
  mutate(Genus = case_when(nSpecies>1~Genus[which.max(X3)],
                           TRUE~Genus)) %>% 
  slice_head() %>% 
  ungroup()

####Look closely at this table to make sure the hits make sense! Maybe before the slice_head step

#now get full taxonomy
lineage_list <- lapply(tax_df$taxid, function(tid) {
  lineage <- entrez_fetch(db="taxonomy", id=tid, rettype="xml")
  # Parse XML to get full hierarchy (Kingdom → Species)
})

parse_lineage <- function(xml_string) {
  doc <- read_xml(xml_string)
  
  # Extract the TaxId for the query
  taxid <- xml_text(xml_find_first(doc, "//Taxon/TaxId"))
  
  # Extract LineageEx Taxon nodes
  lineage_nodes <- xml_find_all(doc, "//LineageEx/Taxon")
  
  # Extract ranks and scientific names
  ranks <- xml_text(xml_find_all(lineage_nodes, "Rank"))
  names <- xml_text(xml_find_all(lineage_nodes, "ScientificName"))
  
  # Combine into a named vector
  lineage <- setNames(names, ranks)
  
  # Add the species itself
  Species <- xml_text(xml_find_first(doc, "//Taxon/ScientificName"))
  lineage["Species"] <- Species
  
  # Return as a data.frame row
  df <- as.data.frame(as.list(lineage), stringsAsFactors=FALSE)
  df$taxid <- taxid
  return(df)
}

taxonomy_df <- lapply(lineage_list, parse_lineage)

# Combine into a single data frame
taxonomy_df <- bind_rows(taxonomy_df)

# Optional: reorder columns: Kingdom → Species
taxonomy_df <- taxonomy_df %>%
  select(kingdom, phylum, class, order, family, genus, Species, taxid) %>% 
  mutate(taxid = as.integer(taxid)) %>% 
  distinct()
colnames(taxonomy_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "taxid")
  
new_taxa <- (unassigned_taxa %>% select(-Genus, -Species)) %>% 
  left_join(taxonomy_df, by = "taxid") %>% 
  select(X1, Kingdom, Phylum, Class, Order, Family, Genus, Species,) %>% 
  column_to_rownames("X1")

compiled_taxa <- as.data.frame(compiled_taxa) %>% 
  filter(!(is.na(Order))) %>% 
  #rename_with(tolower) %>%
  bind_rows(new_taxa) %>% 
  as.matrix()

# Copies from Genus to Species if there is taxa in Genus but NA in species
compiled_taxa <- compiled_taxa %>%
  as.data.frame() %>%
  mutate(Species = ifelse(is.na(Genus) | Genus == "", 
                          Species, 
                          ifelse(is.na(Species) | Species == "", 
                                 Genus, 
                                 Species))) %>%
  as.matrix()


