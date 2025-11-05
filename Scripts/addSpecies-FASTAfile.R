# Code for converting FASTA into a format to use addSpecies()

library(Biostrings)

# 12S
fasta <- readDNAStringSet("DADA2/Ref-DB/MURI_MFU_07_2025.fasta")
headers <- names(fasta)

# Split headers by ";" and extract the last two entries as genus and species
split_headers <- strsplit(headers, ";")

# Get genus (second-to-last) and species (last)
gs_names <- sapply(split_headers, function(x) {
  n <- length(x)
  if (n >= 2) {
    # Remove any extra whitespace
    genus <- trimws(x[n-1])
    species <- trimws(x[n])
    paste(genus, species)
  } else {
    NA  # Mark as NA if not enough taxonomy fields
  }
})

# Filter out incomplete entries (where gs_names is NA)
valid <- !is.na(gs_names)
fasta_gs <- fasta[valid]
gs_names <- gs_names[valid]

# Update headers
names(fasta_gs) <- gs_names

# Makes sure header is formatted correctly
system("sed -E 's/^>([A-Za-z]+) \\1 />\\1 /' DADA2/Ref-DB/ADDSPECIES_MURI_MFU_07_2025.fasta > DADA2/Ref-DB/ADDSPECIES_FIXED.fasta")

# Remove duplicate genus name if present, e.g. "Galeorhinus Galeorhinus galeus" -> "Galeorhinus galeus"
fixed_names <- gsub("^([A-Za-z]+) \\1 ", "\\1 ", names(fasta_gs))

# Apply fixed names to fasta object
names(fasta_gs) <- fixed_names

# Write the fixed fasta for addSpecies
writeXStringSet(fasta_gs, "DADA2/Ref-DB/ADDSPECIES_FIXED.fasta")

# Write to new FASTA
writeXStringSet(fasta_gs, "DADA2/Ref-DB/ADDSPECIES_MURI_MFU_07_2025.fasta")

