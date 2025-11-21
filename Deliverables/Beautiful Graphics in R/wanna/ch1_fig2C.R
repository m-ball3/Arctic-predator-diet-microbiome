library(tidyverse)

PAi <- data.frame(read.csv(file="./Deliverables/Beautiful Graphics in R/wanna/bugs_adj_melt_edit.csv"))
PAi2 <- subset(PAi, PAi$variable!="Coleomegilla")

PAi_sum <- PAi2 %>% 
  group_by(week, year, order, variable) %>% 
  summarise(values = sum(value))

PAi_sum$variable = ordered(PAi_sum$variable,levels=c(
                                         "Other_Coleoptera","Carabidae","Staphylinidae",
                                         "Heteroceridae","Hydrophilidae","Elateridae",
                                         "Phyllophaga","Harmonia","Scarabaeidae",
                                         "Culicidae_and_Chironomidae","Other_Diptera","Tachinidae",
                                         "Muscidae","Sarcophagidae","Syrphidae",
                                         "Cicadellidae","Corixidae","Other_Hemiptera","Miridae",
                                         "Formicidae","Braconidae","Ichneumonidae","Other_Hymenoptera",
                                         "micromoths","Other_Lepidoptera","Noctuidae",
                                         "Arctiidae","Lasiocampidae","Geometridae",
                                         "Sphingidae",
                                         "Araneae", 
                                         "Trichoptera",
                                         "Parasitiformes","Other_Hexapoda","Orthoptera","Plecoptera",
                                         "Opiliones","Neuroptera"))

bug_cols3<-c(
             #coleoptera: 9 shades
            "#003029", "#005045","#007061","#00907C","#00A08A", "#33B3A1", "#66c6b9","#99d9d0","#ccece8",
             #diptera: 6 shades
             "#916800","#c28a00","#F2AD00", "#f5bd33", "#f7ce66","#fade99",
             #hemiptera: 4 shades
             "#2e5e6b","#408496","#5BBCD6","#8cd0e2",
             #hymenoptera: 4 shades
             "#8e7a68","#bda28b","#ECCBAE", "#f2dbc6",
             #leps: 7 shades
             "#01202e","#02364d","#034c6c","#046C9A", "#3689ae", "#68a7c2","#9bc4d7",
             #Araneae
            "#FF0000",
            #trichoptera
             "#D69C4E",
            
             #others: 6 shades
             "#445859","#678585","#89b1b2","#ABDDDE","#d5eeef","#eef8f8")


p2 <- ggplot(PAi_sum, aes(x = week, y = values, fill = variable, na.rm = TRUE))+
  geom_bar(stat = "identity", position ="fill") + facet_wrap(~year + order, ncol=6)
p22 <- p2 + scale_fill_manual(values=bug_cols3)+
  theme_minimal() + theme(legend.position = "bottom",
                        legend.title = element_blank()) 
p22


