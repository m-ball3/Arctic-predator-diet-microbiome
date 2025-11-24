# tries to make the monkey plot
library(ggplot2)
library(dplyr)
library(tidyr)

combined_df <- read.csv("metadata/workable.df.csv")

df_12s <- combined_df %>%                               
  group_by(Predator, Month, Order, Year, Location) %>%                       # group by Predator, Month, Order
  summarise(
    total_asv = sum(Abundance, na.rm = TRUE),                # sum Abundance per cell
    n = n()
  ) %>%
  ungroup()

heatmap <- ggplot(df_12s, aes(x = Location, y = Order, fill = total_asv)) +
  geom_tile(color = "white") +
  facet_wrap(~ Predator, ncol = 1) +
  scale_fill_gradientn(
    colors = c("#ffffd9", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404"),
    breaks = c(0, 1, 5, 10, 15, 20, 25),
    labels = c( "0", "1-5", "6-10", "11-15", "16-20", "21-25", "26+")
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position="bottom")+
  labs(
    x = NULL, y = NULL, fill = "Total ASV Abundance",
    title = "Monthly ASV Abundance by Order and Predator (12S)"
  )
 
heatmap


# 12s only
library(ggplot2)
library(dplyr)
library(tidyr)

df_12s <- combined_df %>%
  filter(Marker == "12S") %>%
  group_by(Predator, Location, Order) %>%
  summarise(
    total_asv = sum(Abundance, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

heatmap.12 <- ggplot(df_12s, aes(x = Location, y = Order, fill = total_asv)) +
  geom_tile(color = "white") +
  facet_wrap(~ Predator, ncol = 1) +
  scale_fill_gradientn(
    colors = c("#ffffd9", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404"),
    breaks = c(0, 1, 5, 10, 15, 20, 25),
    labels = c("0", "1-5", "6-10", "11-15", "16-20", "21-25", "26+")
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom") +
  labs(
    x = NULL, y = NULL, fill = "Total ASV Abundance",
    title = "ASV Abundance by Order & Location (12S)"
  )

heatmap.12















map.nm <- # Combine plots and center legend
  final_plot <- (final_map + nmds.12) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

map.12s.w.faucet.and.nmds <- (final_map / nmds.12.shapes) | faucet.12s

ggsave(
  filename = "./Deliverables/Beautiful Graphics in R/map.12s.w.faucet.and.nmds.png",
  plot = map.12s.w.faucet.and.nmds ,
  width = 20,
  height = 20,
  dpi = 350,
  units = "in",
  bg = "white"
)





























# location NMDS
# Load phyloseq objects
ps.12s.raw <- readRDS("ps.12s.raw")
ps.16s.raw <- readRDS("ps.16s.raw")

# Transform counts for both datasets
ps.prop.12s <- transform_sample_counts(ps.12s.raw, function(otu) otu / sum(otu))
sample_data(ps.prop.12s) <- sample_data(ps.12s.raw)

ps.prop.16s <- transform_sample_counts(ps.16s.raw, function(otu) otu / sum(otu))
sample_data(ps.prop.16s) <- sample_data(ps.16s.raw)
otu_tab <- as.data.frame(otu_table(ps.prop.16s))
ps.prop.16s <- subset_samples(ps.prop.16s, sample_sums(ps.prop.16s) > 0)
ps.prop.16s <- prune_taxa(taxa_sums(ps.prop.16s) > 0, ps.prop.16s)


# Ordinate both objects (you can use same distance for both, e.g. Bray-Curtis)
ord.nmds.12s <- ordinate(ps.prop.12s, method = "NMDS", distance = "bray")
ord.nmds.16s <- ordinate(ps.prop.16s, method = "NMDS", distance = "bray")

# 12S scores/metadata
nmds_12s <- as.data.frame(ord.nmds.12s$points)
nmds_12s$SampleID <- rownames(nmds_12s)
meta_12s <- as.data.frame(sample_data(ps.prop.12s))
meta_12s$SampleID <- rownames(meta_12s)

plot_12s <- left_join(nmds_12s, meta_12s, by = "SampleID") %>%
  mutate(Marker = "12S")

# 16S scores/metadata
nmds_16s <- as.data.frame(ord.nmds.16s$points)
nmds_16s$SampleID <- rownames(nmds_16s)
meta_16s <- as.data.frame(sample_data(ps.prop.16s))
meta_16s$SampleID <- rownames(meta_16s)

plot_16s <- left_join(nmds_16s, meta_16s, by = "SampleID") %>%
  mutate(Marker = "16S")

# Combine datasets
plot_both <- bind_rows(plot_12s, plot_16s)

ggplot(plot_both, aes(x = MDS1, y = MDS2, color = Predator, shape = Location)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Bray NMDS: 12S and 16S",
    x = "NMDS1",
    y = "NMDS2",
    color = "Predator", shape = "Location"
  ) +
  theme_minimal(base_size = 14) +
  facet_wrap(~ Marker)   # If you want panels for each marker




 (final_map + nmdsloc) / heatmap
