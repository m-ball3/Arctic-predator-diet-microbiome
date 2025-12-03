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

# First, check how many unique values there are:
unique(combined_df_with_coords$Location)

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


# Combine plots and center legend
final_plot <- (div.12 + nmds.12) +
  plot_layout(guides = "collect") & theme(legend.position = "none")



final_plot / faucet.12s

ggsave(
  filename = "./Deliverables/Beautiful Graphics in R/final_patchwork_12s-revised.png",
  plot = final_plot / faucet.12s,
  width = 18,
  height = 16,
  dpi = 350,
  units = "in",
  bg = "white"
)


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop.16s <- transform_sample_counts(ps.16s, function(otu) otu/sum(otu))

all(sample_names(ps.prop.16s) %in% rownames(sample_data(ps.prop.16s)))
table(is.na(sample_data(ps.prop.16s)$Predator))

ord_plot_data <- plot_ordination(ps.prop.16s, ord.nmds.bray, color="Predator", title="Bray NMDS")$data
table(is.na(ord_plot_data$Predator))
ord.nmds.bray <- ordinate(ps.prop.16s, method="NMDS", distance="bray")
plot_ordination(ps.prop.16s, ord.nmds.bray, color="Predator", title="Bray NMDS")



library(ggplot2)
ggplot(ord_plot_data[!is.na(ord_plot_data$Predator), ], aes(x=NMDS1, y=NMDS2, color=Predator)) +
  geom_point(size=3, alpha=0.9) +
  labs(title="Bray NMDS (Samples only)") +
  theme_minimal(base_size = 14)
