# Box plot
bp <- ggplot(iris, aes(Species, Sepal.Length)) + 
  geom_boxplot(aes(fill = Species)) +
  theme_minimal() +
  theme(legend.position = "top")
bp