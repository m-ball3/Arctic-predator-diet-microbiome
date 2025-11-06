library(Biostrings)

# Step 1: Read original fasta
fasta <- readDNAStringSet("DADA2/Ref-DB/MURI_MFU_07_2025.fasta")
headers <- names(fasta)

# Extract genus and species from semicolon-delimited headers
split_headers <- strsplit(headers, ";")

# Create Genus_species names
gs_names <- sapply(split_headers, function(x) {
  n <- length(x)
  if (n >= 2) {
    genus <- trimws(x[n-1])
    species <- trimws(x[n])
    paste(genus, species)
  } else {
    NA
  }
})

valid <- !is.na(gs_names)
fasta_gs <- fasta[valid]
gs_names <- gs_names[valid]

# Update headers to Genus_species
names(fasta_gs) <- gs_names

# Step 2: Write intermediate fasta with Genus_species headers
intermediate_file <- "DADA2/Ref-DB/ADDSPECIES_INTERMEDIATE.fasta"
writeXStringSet(fasta_gs, intermediate_file)

# Step 3: Rename sequence IDs to sq1, sq2, etc., preserving Genus_species
lines <- readLines(intermediate_file)
count <- 1
for (i in seq_along(lines)) {
  if (startsWith(lines[i], ">")) {
    parts <- strsplit(lines[i], "\\s+")[[1]]
    if (length(parts) >= 2) {
      # Replace sequence ID with sq{count}, keep species name after that
      lines[i] <- paste0(">sq", count, " ", paste(parts[-1], collapse = " "))
    } else {
      lines[i] <- paste0(">sq", count)
    }
    count <- count + 1
  }
}

# Step 4: Write final fixed fasta
output_file <- "DADA2/Ref-DB/ADDSPECIES_FIXED_FINAL.fasta"
writeLines(lines, output_file)


# Example usage:
# rename_fasta_headers("input.fasta", "output_renamed.fasta")


