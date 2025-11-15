install.packages("writexl")
library(writexl)

ls()
class("taxa")

df <- as.data.frame(merged.taxa)
df_with_rownames <- cbind(rownames = rownames(df), df)
writexl::write_xlsx(df_with_rownames, "WADE003-arcticpred_dada2_QAQC_16SP1+2.Rdata.xlsx")


write_xlsx(df, "mWADE003-arcticpred-16SP2-taxa.xlsx")

getwd()
