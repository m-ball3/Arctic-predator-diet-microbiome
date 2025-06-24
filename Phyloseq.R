# ------------------------------------------------------------------
# FROM DADA2 TO PHYLOSEQ
# ------------------------------------------------------------------


## Sets Working Directory
setwd("C:/Users/MBall/OneDrive/Documents/UW-DOCS/WADE lab/Arctic Predator/DADA2/DADA2 Outputs")

## Libraries

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

## Constructs a dataframe

# creates a character vector of sample IDs
samples.out <- rownames(seqtab.nochim)

# extracts the specific sample ID
# (splits the character at specified ""; extracts the first part)
ID <- sapply(strsplit(samples.out, "003-"), `[`, 2)

# Other options, specific to tutorial
Predator <- substr(ADFG.metadata,1,1)
ADFG.ID <- substr(ADFG.metadata,2,999)
# day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))

# Creates data frame
## HOW IS IT GOING TO KNOW WHICH WADE ID = ADFG ID
samdf <- data.frame(WADE.ID=ID, Predator=Predator, ADFG.ID=ADFG.ID, )

# Get row names from source object
source_rownames <- rownames(seqtab.nochim)
maybe <- data.frame(row.names = source_rownames, taxa=)



# Classifies samples as early or late - probably won't use
# samdf$When <- "Early"
# samdf$When[samdf$Day>100] <- "Late"

# Assigns rownamesaccording to samples.out
rownames(samdf) <- samples.out

# Creates phyloseq object from DADA2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample


# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object,
# and then rename our taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


plot_richness(ps, x="WADE_ID", measures=c("Shannon", "Simpson"), color="When")
