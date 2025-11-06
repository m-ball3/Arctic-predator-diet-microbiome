<<<<<<< HEAD
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
=======
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
>>>>>>> f43aa860c3a2f58828b083a697e6d9bc4fb527d9
    genus <- trimws(x[n-1])
    species <- trimws(x[n])
    paste(genus, species)
  } else {
<<<<<<< HEAD
    NA
  }
})

=======
    NA  # Mark as NA if not enough taxonomy fields
  }
})

# Filter out incomplete entries (where gs_names is NA)
>>>>>>> f43aa860c3a2f58828b083a697e6d9bc4fb527d9
valid <- !is.na(gs_names)
fasta_gs <- fasta[valid]
gs_names <- gs_names[valid]

<<<<<<< HEAD
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

=======
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
>>>>>>> f43aa860c3a2f58828b083a697e6d9bc4fb527d9

