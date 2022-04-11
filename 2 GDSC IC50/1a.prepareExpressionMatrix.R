#' Purpose:
#' To clean and filter GDSC RNA expression data, create "Cell_line_RMA_proc_basalExp.RData"

# Load normalized cell line expression matrix 
exp_matrix <- read.delim(file ="Cell_line_RMA_proc_basalExp.txt", header = T, sep = "\t")
dim(exp_matrix) #dimensions of the matrix = 17419  rows (genes) x 1018 columns (cell)

# Back up expression matrix in an RData file (takes less time to read in)
exp_matrix -> Cell_line_RMA_proc_basalExp
save(Cell_line_RMA_proc_basalExp, file = "Cell_line_RMA_proc_basalExp.RData")
head(exp_matrix[,1:10]) -> x

# Look for duplicated values in gene names
# length(which(exp_matrix$GENE_SYMBOLS == ""))
# y <- as.character(exp_matrix$GENE_title[(which(exp_matrix$GENE_SYMBOLS == ""))])
# y[y!=""]
# z <- exp_matrix$GENE_SYMBOLS
# z <- z[z!=""]
# any(duplicated(z)) #FALSE - no duplicates if there are no gene symbols

# Remove rows with no gene symbols
exp_matrix <- exp_matrix[exp_matrix$GENE_SYMBOLS != "", ]

# Rename heading to cosmic identifiers
# Old names: colnames(exp_matrix)[1:10]
# [1] "DATA.906826"  "DATA.687983"  "DATA.910927"  "DATA.1240138" "DATA.1240139" "DATA.906792"  "DATA.910688"  "DATA.1240135"
colnames(exp_matrix) <- gsub(".*\\.", "", colnames(exp_matrix))
# New names: colnames(exp_matrix)[1:10]
# [1] "906826"  "687983"  "910927"  "1240138" "1240139" "906792"  "910688"  "1240135" "1290812" "907045"

# rename rows as gene names (i.e. ensembl ids) and remove this column from the matrix
rownames(exp_matrix) <- exp_matrix$GENE_SYMBOLS
exp_matrix$GENE_SYMBOLS <- NULL
exp_matrix$GENE_title <- NULL

# check whether NAs exist in matrix using any() 
dim(exp_matrix) #= 17419 rows x 1018 columns
any(is.na(exp_matrix)) # FALSE

# Save clean data
save(exp_matrix, file = "Cell_line_RMA_proc_basalExp_clean.RData")


