library(phyloseq); packageVersion("phyloseq")


## WOULD RATHER HAVE COMBINED AND LOOK AT BY MARKER AND BY PREDATOR???
# Loads in phyloseq objects
ps.12s <- readRDS("ps.12s")
ps.16s <- readRDS("ps.16s")

#alpha richness
plot_richness(ps.12s, x="Year", measures=c("Shannon", "Simpson"), color="Predator")
plot_richness(ps.16s, x="Year", measures=c("Shannon", "Simpson"), color="Predator")


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop.12s <- transform_sample_counts(ps.12s, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop.12s, method="NMDS", distance="bray")
plot_ordination(ps.prop.12s, ord.nmds.bray, color="Predator", title="Bray NMDS")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop.16s <- transform_sample_counts(ps.16s, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop.16s, method="NMDS", distance="bray")
plot_ordination(ps.prop.16s, ord.nmds.bray, color="Predator", title="Bray NMDS")
