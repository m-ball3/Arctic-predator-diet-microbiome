# ------------------------------------------------------------------
# CREATES A TABLE OF UNIQUE SPECIES
# TO BUILD ARCTIC SPECIES CATALOG OF CO1 DATASET
# ------------------------------------------------------------------
library(seqinr)
getwd()

# Creates a list from the .fasta file
fasta_data <- read.fasta("DADA2/Ref-DB/MURI_MFU_07_2025.fasta", 
                         seqtype = "DNA", 
                         as.string = TRUE, 
                         whole.header = TRUE)


# Extracts taxonomic information from the list 
headers <- names(fasta_data)

# Keeps only unique headers and puts them in a data frame
unique<- unique(headers)

split_headers <- strsplit(unique, ";")
taxa_df <- do.call(rbind, lapply(split_headers, function(x) {
  length(x) <- 7
  x
}))
taxa_df <- as.data.frame(taxa_df, stringsAsFactors = FALSE)
colnames(taxa_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
