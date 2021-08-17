# Converts data in excel to RData list
library(openxlsx)

## Functions ------------------------------------------------
# This function reads a list from an excel file
loadList <- function(file){
  library(openxlsx)
  file_df = read.xlsx(file)
  list = lapply(file_df, function(x) x[!is.na(x)])
  return(list)
}

# This function writes a list (named) with vectors of varying sizes to Excel
# It will make a worksheet in the workbook (wb) specified with the name of "current_sheet"
write_list_to_Excel <- function(current_sheet, my_list, wb){
  # # call required package to write to Excel
  require(openxlsx)
  # add a worksheet to the workbook with the name specified in the argument 
  addWorksheet(wb, current_sheet)
  col_num <- 1
  # In a loop, make a column for the signatures' ensembl ids
  for(i in 1:length(my_list)){
    # rename columns to specify cell type
    writeData(wb, sheet = current_sheet, x = names(my_list[i]), startCol = col_num, startRow = 1)
    # write data to worksheet
    writeData(wb, sheet = current_sheet, x = my_list[[i]][1], startCol = col_num, startRow = 2)
    # rename columns to specify cell type
    writeData(wb, sheet = current_sheet, x = names(my_list[i]), startCol = col_num+1, startRow = 1)
    # write data to worksheet
    writeData(wb, sheet = current_sheet, x = my_list[[i]][2], startCol = col_num+1, startRow = 2)
    # update counter
    col_num <- col_num + 2
  }
  return(wb)
}



## MAIN - Run this ------------------------------------------------
# Load signatures from xlsx file
signatures <- loadList(file = "HK_lineage_signatures.xlsx")

# Save as RData file
RData_filename <- "HKs_m_to_hm_5fold_signatures.RData"
save(signatures, file = RData_filename)
