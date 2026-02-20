# ------------------------------------------------------------------
# TRIMLEFT: 30 VS 20 TAXONOMIC ASSIGNMENT COMPARISON
# ------------------------------------------------------------------

# Loads libraries
library(ggplot2)
library(tidyverse)
library(patchwork)

getwd()
trim30 <- readxl::read_excel("WADE003-arcticpred_dada2_QAQC_12SP1_output-130trunc-ADFGnotes-50BOOTTRUE.xlsx")

plot.30 <- ggplot(trim30, aes(x = boot.Species)) +
  geom_histogram()+
  theme_minimal()+
  scale_y_continuous(expand = c(0,0))+
  labs(title = "TrimLeft = 30")+
  theme(axis.text.x = element_blank())+ 
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0, 115))
print(plot.30)

trim20 <- readxl::read_excel("WADE003-arcticpred_dada2_QAQC_12SP1_output-130trunc-ADFGnotes-50BOOTTRUE-leftrim20.xlsx")

plot.20 <- ggplot(trim20, aes(x = boot.Species)) +
  geom_histogram()+
  theme_minimal()+ 
  scale_y_continuous(expand = c(0,0))+
  labs(title = "TrimLeft = 20") + 
  labs(x = "Species Assignment Bootstrap Value")+
  scale_y_continuous(limits = c(0, 115))
print(plot.20)

plots <- plot.30 / plot.20
plots

ggsave("Deliverables/trim-hist.png", plot = plots, width = 10, height = 5, units = "in", dpi = 300)

