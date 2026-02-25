###Filter MIDORI to only
library(Biostrings)

midori <- readDNAStringSet("G:/My Drive/03 Current Research/06 ADFG diet study/Data analysis/CO1 test/MIDORI2_UNIQ_NUC_GB269_CO1_DADA2.fasta.gz")
length(midori)

primerF <- DNAString("GGWACWGGWTGAACWGTWTAYCCYCC")
primerR <- DNAString(chartr("I", "N", "TAIACYTCIGGRTGICCRAARAAYCA"))

# Forward primer matches
fwd_hits <- vmatchPattern(primerF, midori, fixed=FALSE)

# Reverse primer matches (remember we want reverse complement)
rev_hits <- vmatchPattern(reverseComplement(primerR), midori, fixed=FALSE)

# chunk_size <- 100000
# chunks <- split(midori, ceiling(seq_along(midori)/chunk_size))

# keep only sequences with both primers
keep_idx <- which(sapply(fwd_hits, length) > 0 & sapply(rev_hits, length) > 0)
midori_trimmed <- midori[keep_idx]

length(midori_trimmed)

#Trim to primers
trimmed_seqs <- DNAStringSet(
  mapply(function(seq, fwd_match, rev_match) {
    # Skip if no matches
    if (length(fwd_match) == 0 | length(rev_match) == 0) return(NULL)
    
    start_pos <- start(fwd_match)[1]  # first forward match
    end_pos <- end(rev_match)[length(rev_match)]  # last reverse match
    
    # Make sure coordinates are valid
    if (start_pos > end_pos) {
      temp <- start_pos
      start_pos <- end_pos
      end_pos <- temp
    }
    
    subseq(seq, start=start_pos, end=end_pos)
  },
  seq=midori_trimmed,
  fwd_match=fwd_hits[keep_idx],
  rev_match=rev_hits[keep_idx])
)

# Remove NULL entries
trimmed_seqs <- trimmed_seqs[!sapply(trimmed_seqs, is.null)]

#write new fasta
writeXStringSet(trimmed_seqs, "G:/My Drive/03 Current Research/06 ADFG diet study/Data analysis/CO1 test/Midori_Leray_trimmed.fasta")
