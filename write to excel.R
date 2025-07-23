install.packages("writexl")
library(writexl)

ls()
class("taxa")

df <- as.data.frame(taxa)
df_with_rownames <- cbind(rownames = rownames(df), df)
writexl::write_xlsx(df_with_rownames, "WADE003-arcticpred-12SP1-taxa-130trunc2.xlsx")


write_xlsx(df, "mWADE003-arcticpred-16SP2-taxa.xlsx")

getwd()
