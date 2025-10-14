install.packages("writexl")
library(writexl)

ls()
class("taxa")

df <- as.data.frame(taxa)
df_with_rownames <- cbind(rownames = rownames(df), df)
writexl::write_xlsx(df_with_rownames, "WADE003-arcticpred_dada2_QAQC_12SP1_output-130trunc-ADFGnotes-50BOOTTRUE-leftrim20.xlsx")


write_xlsx(df, "mWADE003-arcticpred-16SP2-taxa.xlsx")

getwd()
