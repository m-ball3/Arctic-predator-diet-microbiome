# Final plot 1: 

# FAUCET
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

# OTHERS
library(phyloseq)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Load phyloseq objects
ps.12s.raw <- readRDS("ps.12s.raw")
ps.16s.raw <- readRDS("ps.16s.raw")

# Extract richness data
pr12s_data <- plot_richness(ps.12s.raw, x = "Year", measures = c("Shannon", "Simpson"), color = "Predator")$data
pr16s_data <- plot_richness(ps.16s.raw, x = "Year", measures = c("Shannon", "Simpson"), color = "Predator")$data

pr12s_data$Year <- factor(pr12s_data$Year, levels = sort(unique(pr12s_data$Year)))
pr16s_data$Year <- factor(pr16s_data$Year, levels = sort(unique(pr16s_data$Year)))
dd <- position_dodge(width = 0.6)

# Find y-limits for unified axis
ylim_shannon <- range(c(
  pr12s_data$value[pr12s_data$variable == "Shannon"],
  pr16s_data$value[pr16s_data$variable == "Shannon"]
), na.rm = TRUE)

dd_box   <- position_dodge2(width = 0.7, preserve = "single")
dd_point <- position_jitterdodge(jitter.width = 0.01,
                                 dodge.width  = 0.45)

div.12 <- ggplot(pr12s_data, aes(x = Year, y = value, color = Predator)) +
  geom_boxplot(aes(group = interaction(Year, Predator)),
               position = dd_box,
               fill = NA, outlier.shape = NA,
               width = 0.45, color = "grey40") +
  geom_jitter(position = dd_point, size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  labs(title = "Alpha Diversity", y = "Diversity Index", x = "Year") +
  scale_y_continuous(limits = ylim_shannon, name = "Diversity Index") +
  theme_minimal(base_size = 18) +  # global bump
  theme(
    legend.position = "none",
    axis.text.x  = element_text(hjust = 1, angle = 45, size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 22),
    strip.text   = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

div.12


# Create 12S NMDS plot (legend hidden)
ps.prop.12s <- transform_sample_counts(ps.12s.raw, function(otu) otu / sum(otu))
sample_data(ps.prop.12s) <- sample_data(ps.12s.raw)
ord.nmds.bray <- ordinate(ps.prop.12s, method = "NMDS", distance = "bray")

nmds_scores <- as.data.frame(ord.nmds.bray$points)
nmds_scores$SampleID <- rownames(nmds_scores)
metadata <- as.data.frame(sample_data(ps.prop.12s))
metadata$SampleID <- rownames(metadata)
plot_data <- left_join(nmds_scores, metadata, by = "SampleID")

# Choose distinct shape codes, e.g.
shape_values <- c(21, 22, 23, 24, 25, 3, 4, 8, 7)  # Give each unique Location its own shape

# Assign colors/shapes in your plot:
nmds.12.shapes <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = Predator, shape = Location)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = shape_values) +
  labs(
    title = "Bray NMDS",
    x = "NMDS1", y = "NMDS2",
    color = "Predator", shape = "Location"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "vertical", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


nmds.12.shapes



nmds.12 <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = Predator)) +
  # ellipses first so points are on top
  stat_ellipse(
    aes(fill = Predator),
    geom   = "polygon",
    alpha  = 0.25,
    color  = NA,        # borders off; or set = "black" if you want outlines
    level  = 0.95,      # confidence level
    type   = "t"        # or "norm" depending on your preference
  ) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Bray NMDS",
    x = "NMDS1",
    y = "NMDS2",
    color = "Predator",
    fill  = "Predator"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face  = "bold",
      size  = 22
    ),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16)
  )

nmds.12

library(png)
library(grid)

silhouettes <- list(
  Beluga  = rasterGrob(readPNG("./Deliverables/Beautiful Graphics in R/silhouettes/beluga.png",       native = TRUE), interpolate = TRUE),
  Ringed  = rasterGrob(readPNG("./Deliverables/Beautiful Graphics in R/silhouettes/ringed seal.png", native = TRUE), interpolate = TRUE),
  Bearded = rasterGrob(readPNG("./Deliverables/Beautiful Graphics in R/silhouettes/bearded seal.png",native = TRUE), interpolate = TRUE)
)

# # Example: add one silhouette to the NMDS 12S plot, bottom-right corner
# nmds.12 <- nmds.12 + 
#   annotation_custom(
#     grob = silhouettes$Bearded,
#     xmin = .25, xmax = .75,
#     ymin = 2.25, ymax = 2.75
#   )+
#   annotation_custom(
#     grob = silhouettes$Ringed,
#     xmin = 2.25, xmax = 2.75,
#     ymin = 0.25, ymax = -.25
#   )+
#   annotation_custom(
#     grob = silhouettes$Beluga,
#     xmin = -2, xmax = -1.5,
#     ymin = 0.5, ymax = -.5
#   )
nmds.12

# Turn off legend on faucet plot
faucet.12s <- faucet.12s + theme(legend.position = "none")

# Collect guides only for the top row
final_plot <- (div.12 + nmds.12) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

final_plot / faucet.12s



ggsave(
  filename = "./Deliverables/Beautiful Graphics in R/final_patchwork_12s-revised.png",
  plot = final_plot / faucet.12s,
  width = 22,
  height = 16,
  dpi = 350,
  units = "in",
  bg = "white"
)
