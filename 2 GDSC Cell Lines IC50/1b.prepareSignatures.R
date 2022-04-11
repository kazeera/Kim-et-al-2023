#' Purpose: 
#' To convert signature data from Excel to RData list
library(openxlsx)
library(kazutils)

# This function reads a list from an Excel file. Each column becomes a list element.
loadList <- function(file){
  library(openxlsx)
  # Read Excel file
  file_df = read.xlsx(file)
  # Create list object
  lapply(file_df, function(x) x[!is.na(x)])
}

# Load signatures from xlsx file
signatures <- loadList(file = "HK_lineage_signatures.xlsx")

# Save as RData file
RData_filename <- "HKs_m_to_hm_5fold_signatures.RData"
save(signatures, file = RData_filename)
