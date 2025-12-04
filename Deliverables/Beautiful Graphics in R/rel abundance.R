library(phyloseq)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

library(RColorBrewer)
library(phyloseq)
library(ggplot2)

ps.12s <- readRDS("ps.12s")
# ======== 12S ========
ps.12s <- tax_glom(ps.12s, "Species.y", NArm = FALSE)
# Filters out anything not in Actinopteri
ps.12s <- subset_taxa(ps.12s, Class == "Actinopteri")
nsamples(ps.12s)

# Remove samples with total abundance < 100
ps.12s <- prune_samples(sample_sums(ps.12s) >= 100, ps.12s)
sample_sums(ps.12s)
nsamples(ps.12s)

ps12s.rel <- transform_sample_counts(ps.12s, function(x) { x_rel <- x / sum(x); x_rel[is.nan(x_rel)] <- 0; x_rel })


# --------- Family
families_12s <- sort(unique(as.character(tax_table(ps12s.rel)[,"Family"])))
n_fam_12s <- length(families_12s)
palette_fam_12s <- rep(c(brewer.pal(8,"Accent"), brewer.pal(8,"Set3"), brewer.pal(12, "Set3"), brewer.pal(8, "Set2")), length.out = n_fam_12s)
names(palette_fam_12s) <- families_12s

fam.rel.plot.12s <- plot_bar(ps12s.rel, fill = "Family") +
  scale_fill_manual(values = palette_fam_12s) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), 
        legend.position="bottom")
fam.rel.plot.12s

# --------- Genus
genera_12s <- sort(unique(as.character(tax_table(ps12s.rel)[,"Genus.y"])))
n_gen_12s <- length(genera_12s)
palette_gen_12s <- rep(c(brewer.pal(8,"Pastel1"), brewer.pal(8,"Pastel2"), brewer.pal(12, "Set3"), brewer.pal(8, "Set2")), length.out = n_gen_12s)
names(palette_gen_12s) <- genera_12s

gen.rel.plot.12s <- plot_bar(ps12s.rel, fill = "Genus.y") +
  scale_fill_manual(values = palette_gen_12s) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gen.rel.plot.12s

# --------- Species
species_12s <- sort(unique(as.character(tax_table(ps12s.rel)[,"Species.y"])))
n_sp_12s <- length(species_12s)
palette_sp_12s <- rep(c(brewer.pal(8,"Pastel2"), brewer.pal(8,"Pastel1"), brewer.pal(12, "Set3"), brewer.pal(8, "Set2")),, length.out = n_sp_12s)
names(palette_sp_12s) <- species_12s

sp.rel.plot.12s <- plot_bar(ps12s.rel, fill = "Species.y") +
  scale_fill_manual(values = palette_sp_12s) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
sp.rel.plot.12s

# --------- Facet-wrapped by Predator (Species-level)
faucet.12s <- plot_bar(ps12s.rel, fill = "Family") +
  scale_fill_manual(values = palette_fam_12s) +
  facet_wrap(~ Predator, ncol = 1, scales = "free_x", strip.position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 10, color = "grey60"),
      axis.title.y = element_text(size = 12, color = "grey40"),
      axis.line.y = element_line(color = "grey80"),
      axis.ticks.y = element_line(color = "grey80"),
    legend.position = "bottom", 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
faucet.12s


faucet.12s <- plot_bar(ps12s.rel, fill = "Family") +
  scale_fill_manual(values = palette_fam_12s) +
  facet_wrap(
    ~ Predator,
    nrow = 1,
    scales = "free_x"
  ) +
  theme_minimal(base_size = 18) +  # bump everything
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y  = element_text(size = 14, color = "grey60"),
    axis.title.y = element_text(size = 18, color = "grey40"),
    axis.line.y  = element_line(color = "grey80"),
    axis.ticks.y = element_line(color = "grey80"),
    strip.text   = element_text(size = 16, face = "bold"),  # facet labels
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
faucet.12s




ggsave(
  filename = "./Deliverables/Beautiful Graphics in R/rel.location.png",
  plot = faucet.12s,
  width = 16,
  height = 20,
  dpi = 350,
  units = "in",
  bg = "white"
)


# ======== 16S ========

ps.16s <- tax_glom(ps.16s, "Species.y", NArm = FALSE)
ps16s.rel <- transform_sample_counts(ps.16s, function(x) { x_rel <- x / sum(x); x_rel[is.nan(x_rel)] <- 0; x_rel })

families_16s <- sort(unique(as.character(tax_table(ps16s.rel)[,"Family"])))
n_fam_16s <- length(families_16s)
palette_fam_16s <- rep(c(brewer.pal(12,"Set3"), brewer.pal(8, "Pastel2")), length.out = n_fam_16s)
names(palette_fam_16s) <- families_16s

fam.rel.plot.16s <- plot_bar(ps16s.rel, fill = "Family") +
  scale_fill_manual(values = palette_fam_16s) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
fam.rel.plot.16s

genera_16s <- sort(unique(as.character(tax_table(ps16s.rel)[,"Genus.y"])))
n_gen_16s <- length(genera_16s)
palette_gen_16s <- rep(c(brewer.pal(8,"Pastel1"), brewer.pal(8,"Pastel2")), length.out = n_gen_16s)
names(palette_gen_16s) <- genera_16s

gen.rel.plot.16s <- plot_bar(ps16s.rel, fill = "Genus.y") +
  scale_fill_manual(values = palette_gen_16s) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
gen.rel.plot.16s

species_16s <- sort(unique(as.character(tax_table(ps16s.rel)[,"Species.y"])))
n_sp_16s <- length(species_16s)
palette_sp_16s <- rep(c(brewer.pal(8,"Accent"), brewer.pal(8,"Pastel2")), length.out = n_sp_16s)
names(palette_sp_16s) <- species_16s

sp.rel.plot.16s <- plot_bar(ps16s.rel, fill = "Species.y") +
  scale_fill_manual(values = palette_sp_16s) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
sp.rel.plot.16s

faucet.16s <- plot_bar(ps16s.rel, fill = "Species.y") +
  scale_fill_manual(values = palette_sp_16s) +
  facet_wrap(~ Predator, ncol=1, scales="free_x", strip.position="right") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.5,"lines"),
        axis.title.x = element_text(margin = margin(t=10)))
faucet.16s

