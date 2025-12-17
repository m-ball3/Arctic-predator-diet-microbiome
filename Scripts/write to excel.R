#install.packages("writexl")
library(writexl)

ls()
class("taxa")

df <- as.data.frame(taxa)
df_with_rownames <- cbind(rownames = rownames(df), df)
writexl::write_xlsx(df_with_rownames, "./Deliverables/16S/WADE003-arcticpred_dada2_QAQC_16SP1+2_regions.xlsx")


write_xlsx(df, "/Deliverables/16S/mWADE003-arcticpred-16SP2-taxa.xlsx")

getwd()
