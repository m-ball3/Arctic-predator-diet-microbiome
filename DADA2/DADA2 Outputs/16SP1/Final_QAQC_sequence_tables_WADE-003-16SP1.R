library(phyloseq)
library(Biostrings)
library(tidyverse)
library(vegan)

setwd("G:/My Drive/06 ADFG diet study/Data analysis/QAQC pipeline/16SP1 output")
load("G:/My Drive/06 ADFG diet study/Data analysis/QAQC pipeline/16SP1 output/ADFG_16SP1_dada2output_95boot_102424.Rdata") #from dada2QAQC_laptop.R

# redo taxonomy ------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, "G:/My Drive/06 ADFG diet study/Data analysis/QAQC pipeline/16S_Arctic_predator_reference_database_2024_addcapelinandcod extraction.fasta", 
                       minBoot = 95, multithread=FALSE) 

# check taxonomy in ref database, if needed
# sq <- getSequences("G:/My Drive/06 ADFG diet study/Data analysis/QAQC pipeline/16S_Arctic_predator_reference_database_2024_addcapelinandcod extraction.fasta")
# ids <- names(sq)
# lens <- sapply(strsplit(ids, ";"), length)
# table(lens)
# which(lens != 7)

# sample metadata --------------------------------------------------------------

samdf <- rownames(seqtab.nochim) %>% 
  as.data.frame() %>% 
  select(sampleID = 1) %>% 
  separate(sampleID, into = c("LabID",NA, NA), sep = "_", remove = FALSE) %>% 
  left_join(read.csv("G:/My Drive/06 ADFG diet study/Metadata/ADFG_dDNA_labwork_metadata.csv"), by = c("LabID" = "LabID")) %>% 
  select(sampleID, LabID, Specimen.ID, Shipment) %>% 
  left_join(read.csv("G:/My Drive/06 ADFG diet study/Metadata/ADFG_dDNA_sample_metadata.csv"), by = c("Specimen.ID" = "Specimen.ID", "Shipment" = "Shipment")) %>% 
  column_to_rownames("sampleID")

# create phyloseq object -------------------------------------------------------

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
# update ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

ps.merge <- tax_glom(ps, "Species", NArm = FALSE)

reads.per.samp <- as.data.frame(sample_sums(ps.merge), nm = c("w_mammals"))

hist(sample_sums(ps.merge), main="Read Counts w mammals", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=6)

# remove hosts, merge to species -----------------------------------------------

ps <- subset_taxa(ps, Class!="Mammalia")
ps <- prune_samples(sample_sums(ps)>=20, ps) #remove samples with < 20 reads
ps.sp <- tax_glom(ps, "Species", NArm = FALSE)
ps.sp <- prune_samples(sample_sums(ps.sp)>=2, ps.sp) #remove samples with < 20 reads

hist(sample_sums(ps.sp), main="Read Counts w.o mammals", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=6)

reads.per.samp <- reads.per.samp %>% 
  mutate(wo_mammals = sample_sums(ps.sp)) %>% 
  mutate(read.diff = (w_mammals - wo_mammals))

# convert to proportion --------------------------------------------------------

ps.sp <- transform_sample_counts(ps.sp, function(OTU) OTU/sum(OTU) )

# remove species < 1% in all samples -------------------------------------------
### not used
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.sp, f1, A=1)
ps.sp.major <- prune_taxa(lowcount.filt, ps.sp)

# Quick PERMANOVA -------------------------------------------------------------

allsamps.otu <- as.data.frame(ps.sp@otu_table) 


tax_table <- as.data.frame(ps.sp@tax_table)

#colnames(allsamps.otu) <- tax_table$Species

sex_difference <- adonis2(allsamps.otu ~ Sex, data = samdf)
sex_difference

bar.sp <- plot_bar(ps.sp, fill = "Species", x = "Specimen.ID", title = "Taxonomy by species")
bar.gen <- plot_bar(ps.sp, fill = "Genus", x = "Specimen.ID", title = "Taxonomy by genus")
  
write.csv(allsamps.otu, file = "ADFG_16SP1_Autumn2024_speciesprop_boot95.csv")
write.csv(tax_table, file = "ADFG_16SP1_Autumn2024_taxonomy_boot95.csv")
write.csv(samdf, file = "ADFG_16SP1_Autumn2024_metadata_boot95.csv")

save(allsamps.otu, tax_table, samdf, ps.sp, file = "ADFG_16SP1_Autumn2024_data_boot95.Rdata")

ggsave("ADFG_16SP1_proportional_diet_species.png", bar.sp, height = 6, width = 8, units = "in")
ggsave("ADFG_16SP1_proportional_diet_genus.png", bar.gen, height = 6, width = 8, units = "in")


refseqs <- as.data.frame(refseq(ps.sp))

write.csv(reads.per.samp, "Read_depth_per_sample_16SP1_121824.csv")
