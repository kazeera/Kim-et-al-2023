# rm(list=ls())
# This script subsets the IC50 matrix to cell lines and drugs of interest

##EDIT on 20191128: Hyeyeon gave me a new drug list with drug identifiers (numerical) to match to the GDSC data
library(openxlsx)
tissue_type <- "breast"
# IC50 matrix holds IC50 data for DRUGS (columns) by CELL LINES (rows)
#load IC50 matrix where row names are cell line identifiers and column names are drug identifiers
IC50_mat <- read.xlsx("logIC50_values_cell_Lines_clean.xlsx", colNames = T, rowNames = T)
IC50_mat$Sample.Names <- NULL

# LOOK AT DRUGS FIRST --------------------------------------------
# 1. MAKE DRUGS OF INTEREST IDENTIFIERS --------------------------
# Read in drug annotations
# a) Read in drug list from Excel file and order lowest to highest identifier
drug_list <- read.xlsx("drug_list_for_gdsc.xlsx") #subset of "TableS1F Drug Ann clean.xlsx"
drug_list <- drug_list[order(drug_list$Identifier),]


# 2. SUBSET IC50 MATRIX to drugs of interest --------------------------------------
#subset original data table based on column names that match the cols interest list
subset_table_cols <- function(table, cols_interest){
  all_cols <- colnames(table)
  table <- table[,all_cols %in% cols_interest]
  colnames(table) <- all_cols[all_cols %in% cols_interest]
  return (table)
}

# Do it. subset matrix by matching columns
IC50_subset <- subset_table_cols(IC50_mat, drug_list$Identifier)
IC50_subset <- IC50_subset[,order(as.numeric(colnames(IC50_subset)))]
all(IC50_mat$'152'== IC50_subset$'152', na.rm = T) # TRUE


# LOOK AT CELL LINES ONLY -----------------------------------------------
# Which ones are tissue type of interest cell lines? did this in step 2 already
cell_lines_of_interest <- colnames(ssGseaValues)

# a) Check that all cell lines in IC50 matrix have annotations
all(rownames(IC50_subset) == rownames(IC50_mat)) # TRUE

# b) subset to cell lines of interest by matching rows
#returns a subset of a "table" with only the rows of interest
#subset original data table based on rownames that match the interest list
subset_table_rows <- function(table, rows_interest){
  all_rows <- rownames(table)
  table <- table[ all_rows %in% rows_interest,]
  rownames(table) <-  all_rows[ all_rows %in% rows_interest]
  return(table)
}

# Do it..
IC50_subset <- subset_table_rows(IC50_subset, cell_lines_of_interest)

# IC50 matrix holds IC50 data for DRUGS OF INTEREST (columns) by CELL LINES OF INTEREST (rows)
head(IC50_subset[,1:6])

save(IC50_subset, file = "IC50_breast.cell.lines_HK.drugs.RData")

