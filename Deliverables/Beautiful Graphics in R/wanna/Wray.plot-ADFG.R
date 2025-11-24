library(tidyverse)
library(dplyr)

df <- data.frame(read.csv(file="metadata/workable.df.csv")) 
  # Abundance = relative abundance by marker
#PAi2 <- subset(PAi, PAi$variable!="Coleomegilla")

df <- df %>% filter(!is.na(Species.y))
 df <- df %>% filter(Abundance > 0)

sp.fam <- df %>%
  select(Species.y, Family, Abundance, Predator) %>%  # Select only the needed columns
  distinct(Species.y, .keep_all = TRUE)   # Keep unique Species.y rows, preserving and Family

# Saves combined_df as a .csv for Amy's abundance table
write.csv(sp.fam, "Deliverables/Beautiful Graphics in R/wanna/sp.fam.csv")
df <- df %>%
  mutate(Season = case_when(
    Month %in% c("DEC", "JAN", "FEB") ~ "Winter",
    Month %in% c('MAR', "APR", 'MAY') ~ "Spring",
    Month %in% c('JUN', 'JUL', 'AUG') ~ "Summer",
    Month %in% c('SEP', 'OCT', 'NOV') ~ "Fall",
    TRUE ~ NA_character_
  )
  )

df_sum <- df %>% 
  group_by(Predator, Month, Year, Family, Species.y, Marker, Abundance, Season)

df_sum$Species.y = ordered(df_sum$Species.y,levels=c(
  # #Salmonidae: 9 shades
  # "Coregonus spp.","Oncorhynchus keta","Oncorhynchus gorbuscha",
  # "Oncorhynchus kisutch","Coregonus laurettae","Coregonus nasus",
  # "Oncorhynchus spp.","Salvelinus spp.",
  #Cottidae: 5 shades
  "Gymnocanthus tricuspis",
  "Myoxocephalus spp.","Gymnocanthus spp.","Hemilepidotus spp.","Gymnocanthus pistilliger",
  #Pleuronectidae: 5 shades
  "Glyptocephalus zachirus","Platichthys stellatus",
  "Hippoglossoides elassodon","Pleuronectes quadrituberculatus","Limanda aspera",
  #Stichaeidaea: 4 shades
  "Acantholumpenus mackayi",
  "Leptoclinus maculatus","Stichaeus punctatus",
  "Lumpenus spp.",
  #Gadidae: 4 shades
  "Eleginus gracilis","Gadus chalcogrammus",
  "Boreogadus saida", 
  #Osmeridae: 3 shades
  "Hypomesus olidus","Osmerus mordax", "Mallotus villosus", 
  #Gasterosteidae: 1 shades
  "Pungitius pungitius",
  #Liparidae: 3 shades
  "Liparis spp.", "Liparis microstomus", "Liparis tunicatus",
  #Other: 4 shades
  "Pallasina barbata","Ammodytes hexapterus","Clupea pallasii","Nautichthys pribilovius"))

fish_cols3<-c(
  # #Salmonidae: 8 shades
  # "#003029", "#005045","#007061","#00907C","#00A08A", "#33B3A1", "#66c6b9","#99d9d0",
  #Cottidae: 9 shades
  "#916800","#c28a00","#F2AD00", "#f5bd33", "#f7ce66","#fade99",
  #Pleuronectidae: 5 shades
  "#2e5e6b","#408496","#5BBCD6","#8cd0e2","#b4d9e9", 
  #Stichaeidaea: 5 shades
  "#8e7a68","#bda28b","#ECCBAE", "#f2dbc6",
  #Gadidae: 4 shades
  "#01202e","#02364d","#034c6c",#"#046C9A", 
  #Osmeridae: 3 shades
  "#3689ae", "#68a7c2","#9bc4d7",
  #Gasterosteidae: 3 shades
  "#CC0000",# "#FF0000","#FF6666" ,
  #Liparidae: 3 shades
  "#A9732A" , "#D69C4E","#F1C77D",
  
  #NA: 6 shades
  "#445859","#678585","#89b1b2","#ABDDDE")#"#d5eeef","#eef8f8")


p2 <- ggplot(df_sum, aes(x = Specimen.ID, y = Abundance, fill = Species.y, na.rm = TRUE)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Predator ~ Marker + Family) +  # Facet grid: rows by Predator, cols by Marker+Family
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#labs(x = NULL) +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p22 <- p2 + 
  scale_fill_manual(values = fish_cols3) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 # labs(x = NULL) +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p22




