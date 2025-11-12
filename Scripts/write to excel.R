install.packages("writexl")
library(writexl)

ls()
class("taxa")

df <- as.data.frame(merged.taxa)
df_with_rownames <- cbind(rownames = rownames(df), df)
writexl::write_xlsx(df_with_rownames, "WADE003-arcticpred_dada2_QAQC_16SP2_output-Addspecies2.xlsx")


write_xlsx(df, "mWADE003-arcticpred-16SP2-taxa.xlsx")

getwd()
