  # ------------------------------------------------------------------
  # BLASTS UNASSIGNED TAXA TO GET BETTER RESOLUTION
  # THIS IS THE FOURTH STEP AFTER DADA2
  ## 1_rownames-match_12S.R
  ## 2_replicates_contaminated_12s.R and
  ## 3_taxonomy-by-region_12S.R must be run before this
  # ------------------------------------------------------------------
  
  # ------------------------------------------------------------------
  # Sets up the Environment and Loads in data
  # ------------------------------------------------------------------
  
  library(Biostrings)
  library(tidyverse)
  library(rentrez)
  library(purrr)
  library(progress)
  library(xml2)
  library(dplyr)
  
  getwd()
  load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-COOKINLET.Rdata")
  load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-SBERING.Rdata")
  load("DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-ARCTIC.Rdata")
  
  # Extract unassigned sequences and write FASTA-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  # Find sequences with NA below Genus (or chosen rank)
  cook.unassigned_idx <- which(is.na(taxa.cook[, "Genus"]))
  
  # Extract sequence names
  cook.seqs_unassigned <- colnames(cook.seqtab)[cook.unassigned_idx]
  
  # Convert to DNAStringSet for export
  cook.seqs_unassigned_fasta <- DNAStringSet(cook.seqs_unassigned)
  names(cook.seqs_unassigned_fasta) <- cook.seqs_unassigned  # sequence itself as ID
  
  writeXStringSet(
    cook.seqs_unassigned_fasta,
    "./DADA2/Ref-DB/12S/reBLAST/unassigned_seqs-COOKINLET.fasta"
  )
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  # Find sequences with NA below Genus (or chosen rank)
  sbering.unassigned_idx <- which(is.na(taxa.sbering[, "Genus"]))
  
  # Extract sequence names
  sbering.seqs_unassigned <- colnames(sbering.seqtab)[sbering.unassigned_idx]
  
  # Convert to DNAStringSet for export
  sbering.seqs_unassigned_fasta <- DNAStringSet(sbering.seqs_unassigned)
  names(sbering.seqs_unassigned_fasta) <- sbering.seqs_unassigned  # sequence itself as ID
  
  writeXStringSet(
    sbering.seqs_unassigned_fasta,
    "./DADA2/Ref-DB/12S/reBLAST/unassigned_seqs-SBERING.fasta"
  )
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  # Find sequences with NA below Genus (or chosen rank)
  arctic.unassigned_idx <- which(is.na(taxa.arctic[, "Genus"]))
  
  # Extract sequence names
  arctic.seqs_unassigned <- colnames(arctic.seqtab)[arctic.unassigned_idx]
  
  # Convert to DNAStringSet for export
  arctic.seqs_unassigned_fasta <- DNAStringSet(arctic.seqs_unassigned)
  names(arctic.seqs_unassigned_fasta) <- arctic.seqs_unassigned  # sequence itself as ID
  
  writeXStringSet(
    arctic.seqs_unassigned_fasta,
    "./DADA2/Ref-DB/12S/reBLAST/unassigned_seqs-ARCTIC.fasta"
  )
  
  
  ## PUT THIS OUTPUT INTO BLAST AND SAVE A .CSV
  
  # Read BLAST hits and pick best hit per query-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  cook.unassigned_hits <- readr::read_csv(
    "./DADA2/Ref-DB/12S/reBLAST/COOKINLET-Alignment-HitTable.csv",
    col_names = FALSE
  )
  
  # Unique accession IDs from BLAST (column X2)
  cook.accessions <- unique(cook.unassigned_hits$X2)
  
  # Safe wrapper around entrez_summary
  safe_esummary <- purrr::safely(function(acc) {
    rentrez::entrez_summary(db = "nuccore", id = acc)
  })
  
  # Helper: fetch accession metadata with progress bar
  fetch_accession_metadata <- function(accessions, sleep = 0.34) {
    accessions <- unique(accessions)
    n_acc <- length(accessions)
    
    pb <- progress::progress_bar$new(
      format = "  accession metadata [:bar] :current/:total (:percent) eta: :eta",
      total = n_acc,
      clear = FALSE,
      width = 80
    )
    
    purrr::map_df(accessions, function(acc) {
      Sys.sleep(sleep)   # polite rate limiting
      pb$tick()
      
      res <- safe_esummary(acc)
      
      if (!is.null(res$error)) {
        message("Failed for accession: ", acc, " | error: ", res$error$message)
        return(tibble(
          accession = acc,
          title     = NA_character_,
          taxid     = NA_integer_,
          organism  = NA_character_
        ))
      }
      
      s <- res$result
      tibble(
        accession = acc,
        title     = s$title    %||% NA_character_,
        taxid     = s$taxid    %||% NA_integer_,
        organism  = s$organism %||% NA_character_
      )
    })
  }
  
  
  # # Get accession-level metadata (title, taxid, organism)
  cook.accessions <- unique(cook.unassigned_hits$X2)
  cook.tax_df <- fetch_accession_metadata(cook.accessions)
  
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  sbering.unassigned_hits <- readr::read_csv(
    "./DADA2/Ref-DB/12S/reBLAST/SBERING-Alignment-HitTable.csv",
    col_names = FALSE
  )
  
  # Unique accession IDs from BLAST (column X2)
  sbering.accessions <- unique(sbering.unassigned_hits$X2)
  
  # Get accession-level metadata (title, taxid, organism)
  sbering.accessions <- unique(sbering.unassigned_hits$X2)
  sbering.tax_df <- fetch_accession_metadata(sbering.accessions)
  
  
  #-------------------------------------------------------------------------------------
  # ARCTIC
  #-------------------------------------------------------------------------------------
  
  arctic.unassigned_hits <- readr::read_csv(
    "./DADA2/Ref-DB/12S/reBLAST/ARCTIC-Alignment-HitTable.csv",
    col_names = FALSE
  )
  
  # Unique accession IDs from BLAST (column X2)
  arctic.accessions <- unique(arctic.unassigned_hits$X2)
  
  # Get accession-level metadata (title, taxid, organism)
  arctic.accessions <- unique(arctic.unassigned_hits$X2)
  arctic.tax_df <- fetch_accession_metadata(arctic.accessions)
  
  # Fold the two versions of unassigned_taxa into one-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  
  # Attach organism/taxid, pick best hit per query, derive genus/species label
  cook.unassigned_taxa <- cook.unassigned_hits %>%
    left_join(cook.tax_df, by = c("X2" = "accession")) %>%
    filter(!is.na(organism)) %>%
    group_by(X1) %>%   # X1 = query sequence ID
    mutate(
      organism_clean = stringr::str_extract(organism, "^[A-Z][a-z]+"),  # genus word
      # Choose genus from hit with max bit score (X3), if multiple genera
      Genus = ifelse(
        dplyr::n_distinct(organism_clean, na.rm = TRUE) > 1,
        organism_clean[which.max(as.numeric(X3))],
        organism_clean[1]
      ),
      Species = "sp."
    ) %>%
    slice_head(n = 1) %>%   # keep one row per query (the top hit)
    ungroup() %>%
    select(X1, Genus, Species, organism, taxid, title, everything())
  
  head(cook.unassigned_taxa)
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  
  # Attach organism/taxid, pick best hit per query, derive genus/species label
  sbering.unassigned_taxa <- sbering.unassigned_hits %>%
    left_join(sbering.tax_df, by = c("X2" = "accession")) %>%
    filter(!is.na(organism)) %>%
    group_by(X1) %>%   # X1 = query sequence ID
    mutate(
      organism_clean = stringr::str_extract(organism, "^[A-Z][a-z]+"),  # genus word
      # Choose genus from hit with max bit score (X3), if multiple genera
      Genus = ifelse(
        dplyr::n_distinct(organism_clean, na.rm = TRUE) > 1,
        organism_clean[which.max(as.numeric(X3))],
        organism_clean[1]
      ),
      Species = "sp."
    ) %>%
    slice_head(n = 1) %>%   # keep one row per query (the top hit)
    ungroup() %>%
    select(X1, Genus, Species, organism, taxid, title, everything())
  
  head(sbering.unassigned_taxa)
  
  
  #-------------------------------------------------------------------------------------
  # ARCTIC
  #-------------------------------------------------------------------------------------
  
  
  # Attach organism/taxid, pick best hit per query, derive genus/species label
  arctic.unassigned_taxa <- arctic.unassigned_hits %>%
    left_join(arctic.tax_df, by = c("X2" = "accession")) %>%
    filter(!is.na(organism)) %>%
    group_by(X1) %>%   # X1 = query sequence ID
    mutate(
      organism_clean = stringr::str_extract(organism, "^[A-Z][a-z]+"),  # genus word
      # Choose genus from hit with max bit score (X3), if multiple genera
      Genus = ifelse(
        dplyr::n_distinct(organism_clean, na.rm = TRUE) > 1,
        organism_clean[which.max(as.numeric(X3))],
        organism_clean[1]
      ),
      Species = "sp."
    ) %>%
    slice_head(n = 1) %>%   # keep one row per query (the top hit)
    ungroup() %>%
    select(X1, Genus, Species, organism, taxid, title, everything())
  
  head(arctic.unassigned_taxa)
  
  
  # Batch taxonomy lineage fetch in a helper-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  
  fetch_lineages <- function(taxids, batch_size = 50, sleep = 0.4) {
    taxids <- taxids[!is.na(taxids)]
    n_taxids <- length(taxids)
    n_batches <- ceiling(n_taxids / batch_size)
    
    cat("Fetching lineage for", n_taxids, "taxids in", n_batches, "batches\n")
    
    pb <- progress_bar$new(total = n_taxids)
    lineage_xml <- vector("list", n_taxids)
    
    idx <- 1
    
    for (b in seq_len(n_batches)) {
      start_idx <- (b - 1) * batch_size + 1
      end_idx   <- min(b * batch_size, n_taxids)
      batch_taxids <- taxids[start_idx:end_idx]
      
      batch_xml <- map(batch_taxids, function(tid) {
        Sys.sleep(sleep)  # rate limiting
        pb$tick()
        tryCatch(
          rentrez::entrez_fetch(db = "taxonomy", id = tid, rettype = "xml"),
          error = function(e) NA_character_
        )
      })
      
      lineage_xml[start_idx:end_idx] <- batch_xml
    }
    
    pb$terminate()
    
    tibble(
      taxid   = taxids,
      lineage = lineage_xml
    )
  }
  
  # Use it:
  cook.lineage_df <- fetch_lineages(cook.tax_df$taxid)
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  # Use it:
  sbering.lineage_df <- fetch_lineages(sbering.tax_df$taxid)
  
  #-------------------------------------------------------------------------------------
  # ARCTIC
  #-------------------------------------------------------------------------------------
  
  # Use it:
  arctic.lineage_df <- fetch_lineages(arctic.tax_df$taxid)
  
  
  # Parse lineage XML to ranks-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  # Helper to parse one lineage XML string into a tibble of ranks
  parse_lineage <- function(taxid, xml_txt) {
    if (is.na(xml_txt) || xml_txt == "") {
      return(tibble(
        taxid  = taxid,
        Kingdom = NA_character_,
        Phylum  = NA_character_,
        Class   = NA_character_,
        Order   = NA_character_,
        Family  = NA_character_,
        Genus   = NA_character_,
        Species = NA_character_
      ))
    }
    
    doc <- read_xml(xml_txt)
    
    # NCBI taxonomy XML has Taxon elements with Rank + ScientificName
    nodes <- xml_find_all(doc, ".//Taxon")
    ranks <- xml_text(xml_find_all(nodes, "./Rank"))
    names <- xml_text(xml_find_all(nodes, "./ScientificName"))
    
    df <- tibble(rank = ranks, name = names)
    
    # Extract the ranks we care about
    get_rank <- function(r) {
      df %>% filter(rank == r) %>% slice_tail(n = 1) %>% pull(name) %>% { if (length(.) == 0) NA_character_ else . }
    }
    
    tibble(
      taxid   = taxid,
      Kingdom = get_rank("superkingdom"),  # often "Bacteria", "Metazoa", etc.
      Phylum  = get_rank("phylum"),
      Class   = get_rank("class"),
      Order   = get_rank("order"),
      Family  = get_rank("family"),
      Genus   = get_rank("genus"),
      Species = get_rank("species")
    )
  }
  
  # Apply to all taxids in cook.lineage_df
  cook.lineage_ranks <- map2_df(
    cook.lineage_df$taxid,
    cook.lineage_df$lineage,
    parse_lineage
  )
  
  head(cook.lineage_ranks)
  
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  # Apply to all taxids in cook.lineage_df
  sbering.lineage_ranks <- map2_df(
    sbering.lineage_df$taxid,
    sbering.lineage_df$lineage,
    parse_lineage
  )
  
  head(sbering.lineage_ranks)
  
  
  #-------------------------------------------------------------------------------------
  # ARCTIC
  #-------------------------------------------------------------------------------------
  
  # Apply to all taxids in cook.lineage_df
  arctic.lineage_ranks <- map2_df(
    arctic.lineage_df$taxid,
    arctic.lineage_df$lineage,
    parse_lineage
  )
  
  head(arctic.lineage_ranks)
  
  # Attach lineage ranks to unassigned ASVs-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  # Join NCBI lineage ranks to BLAST-unassigned taxa
  cook.unassigned_taxa_full <- cook.unassigned_taxa %>%
    left_join(cook.lineage_ranks, by = "taxid")
  
  head(cook.unassigned_taxa_full)
  
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  # Join NCBI lineage ranks to BLAST-unassigned taxa
  sbering.unassigned_taxa_full <- sbering.unassigned_taxa %>%
    left_join(sbering.lineage_ranks, by = "taxid")
  
  head(sbering.unassigned_taxa_full)
  
  #-------------------------------------------------------------------------------------
  # ARCTIC
  #-------------------------------------------------------------------------------------
  
  # Join NCBI lineage ranks to BLAST-unassigned taxa
  arctic.unassigned_taxa_full <- arctic.unassigned_taxa %>%
    left_join(arctic.lineage_ranks, by = "taxid")
  
  head(arctic.unassigned_taxa_full)
  
  # Merge back into taxa-----------------------------------------
  
  #-------------------------------------------------------------------------------------
  # COOK INLET
  #-------------------------------------------------------------------------------------
  
  
  # Turn taxa.cook into a tibble with an ASV ID column
  taxa_cook_df <- taxa.cook %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ASV")
  
  # Merge NCBI lineage for unassigned ASVs
  taxa_cook_updated <- taxa_cook_df %>%
    left_join(
      cook.unassigned_taxa_full %>%
        select(ASV = X1, Phylum, Class, Order, Family, Genus.y, Species.y),
      by = "ASV",
      suffix = c(".old", ".new")
    ) %>%
    mutate(
      # For each rank, prefer new NCBI lineage where old was NA
      # Kingdom = if_else(is.na(Kingdom.old) & !is.na(Kingdom.new), Kingdom.new, Kingdom.old),
      Phylum  = if_else(is.na(Phylum.old)  & !is.na(Phylum.new),  Phylum.new,  Phylum.old),
      Class   = if_else(is.na(Class.old)   & !is.na(Class.new),   Class.new,   Class.old),
      Order   = if_else(is.na(Order.old)   & !is.na(Order.new),   Order.new,   Order.old),
      Family  = if_else(is.na(Family.old)  & !is.na(Family.new),  Family.new,  Family.old),
      Genus   = if_else(is.na(Genus)   & !is.na(Genus.y),   Genus.y,   Genus),
      Species = if_else(is.na(Species) & !is.na(Species.y), Species.y, Species)
    ) %>%
    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, DB)
  
  
  
  #-------------------------------------------------------------------------------------
  # S BERING SEA
  #-------------------------------------------------------------------------------------
  
  
  # Turn taxa.cook into a tibble with an ASV ID column
  taxa_sbering_df <- taxa.sbering %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ASV")
  
  # Merge NCBI lineage for unassigned ASVs
  taxa_sbering_updated <- taxa_sbering_df %>%
    left_join(
      sbering.unassigned_taxa_full %>%
        select(ASV = X1, Phylum, Class, Order, Family, Genus.y, Species.y),
      by = "ASV",
      suffix = c(".old", ".new")
    ) %>%
    mutate(
      # For each rank, prefer new NCBI lineage where old was NA
      # Kingdom = if_else(is.na(Kingdom.old) & !is.na(Kingdom.new), Kingdom.new, Kingdom.old),
      Phylum  = if_else(is.na(Phylum.old)  & !is.na(Phylum.new),  Phylum.new,  Phylum.old),
      Class   = if_else(is.na(Class.old)   & !is.na(Class.new),   Class.new,   Class.old),
      Order   = if_else(is.na(Order.old)   & !is.na(Order.new),   Order.new,   Order.old),
      Family  = if_else(is.na(Family.old)  & !is.na(Family.new),  Family.new,  Family.old),
      Genus   = if_else(is.na(Genus)   & !is.na(Genus.y),   Genus.y,   Genus),
      Species = if_else(is.na(Species) & !is.na(Species.y), Species.y, Species)
    ) %>%
    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, DB)
  
  #-------------------------------------------------------------------------------------
  # ARCTIC
  #-------------------------------------------------------------------------------------
  
  
  # Turn taxa.cook into a tibble with an ASV ID column
  taxa_arctic_df <- taxa.arctic %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ASV")
  
  # Merge NCBI lineage for unassigned ASVs
  taxa_arctic_updated <- taxa_arctic_df %>%
    left_join(
      arctic.unassigned_taxa_full %>%
        select(ASV = X1, Phylum, Class, Order, Family, Genus.y, Species.y),
      by = "ASV",
      suffix = c(".old", ".new")
    ) %>%
    mutate(
      # For each rank, prefer new NCBI lineage where old was NA
      # Kingdom = if_else(is.na(Kingdom.old) & !is.na(Kingdom.new), Kingdom.new, Kingdom.old),
      Phylum  = if_else(is.na(Phylum.old)  & !is.na(Phylum.new),  Phylum.new,  Phylum.old),
      Class   = if_else(is.na(Class.old)   & !is.na(Class.new),   Class.new,   Class.old),
      Order   = if_else(is.na(Order.old)   & !is.na(Order.new),   Order.new,   Order.old),
      Family  = if_else(is.na(Family.old)  & !is.na(Family.new),  Family.new,  Family.old),
      Genus   = if_else(is.na(Genus)   & !is.na(Genus.y),   Genus.y,   Genus),
      Species = if_else(is.na(Species) & !is.na(Species.y), Species.y, Species)
    ) %>%
    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, DB)
  
  
  
  # save.image(file = "post-BLAST-unassigned.RData")
  save(samdf, cook.seqtab, freq.nochim, track, taxa.cook, taxa_cook_updated, file = "DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-COOKINLET.Rdata")
  save(samdf, sbering.seqtab, freq.nochim, track, taxa.sbering,taxa_sbering_updated, file = "DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-SBERING.Rdata")
  save(samdf, arctic.seqtab, freq.nochim, track, taxa.arctic,taxa_arctic_updated, file = "DADA2/DADA2 Outputs/WADE003-arcticpred_dada2_QAQC_12SP1_output-ARCTIC.Rdata")
  
  
  
  
  # OLD CODE
  
  # cook.tax_df <- map_df(cook.accessions, function(acc) {
  #   res <- safe_esummary(acc)
  #   if (!is.null(res$error)) {
  #     message("Failed for accession: ", acc, " | error: ", res$error)
  #     return(tibble(
  #       accession = acc,
  #       title     = NA_character_,
  #       taxid     = NA_integer_,
  #       organism  = NA_character_
  #     ))
  #   }
  #   s <- res$result
  #   tibble(
  #     accession = acc,
  #     title     = s$title      %||% NA_character_,
  #     taxid     = s$taxid      %||% NA_integer_,
  #     organism  = s$organism   %||% NA_character_
  #   )
  # })
  
  
  
