#install.packages("writexl")

library(writexl)
library(dplyr)


# full tax table with rownames as a column
df_raw <- as.data.frame(taxa.sbering)
df_raw <- df_raw %>%
  tibble::rownames_to_column("ASV")

# Raw: everything
df_raw_sheet <- df_raw

# Filt: Class == Actinopterygii
df_filt_sheet <- df_raw %>%
  filter(Class == "Actinopteri")

# Unassigned: Species is NA or blank
df_unassigned_sheet <- df_filt_sheet %>%
  filter(is.na(Species) | Species == "")

# spp.: Species contains "spp."
df_spp_sheet <- df_filt_sheet %>%
  filter(grepl("spp\\.", Species))

# Write all four to one workbook
write_xlsx(
  list(
    Raw        = df_raw_sheet,
    Filt       = df_filt_sheet,
    Unassigned = df_unassigned_sheet,
    spp        = df_spp_sheet
  ),
  path = "./Deliverables/12S/regions/sbering/WADE003-arcticpred_dada2_QAQC_12SP1_sbering.xlsx"
)


## OLD CODE
# ls()
# class("taxa")
# 
# df <- as.data.frame(taxa.arctic)
# df_with_rownames <- cbind(rownames = rownames(df), df)
# writexl::write_xlsx(df_with_rownames, "./Deliverables/12S/regions/arctic/WADE003-arcticpred_dada2_QAQC_12SP1_arctic.xlsx")
# 

