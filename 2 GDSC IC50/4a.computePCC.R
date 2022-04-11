#' Purpose:
#' To compute PCC (pearson correlation coefficients) of IC50 values and ssGSEA scores
#' 
#' Credit: Pujana lab for conception, Luis for PCC calculations

# Load required R libraries
library(foreach) # to do drug correlations
library(reshape2) # dcast - table reformatting function
library(pheatmap) # heatmap visualization
library(RColorBrewer) # heatmap color
library(scales) # rescale matrix

# Define parameters of analysis (defined in script 0 also, please change for other analyses)
output_folder <- "Human Breast Cell Lines IC50 (FC5, p0.05, human_signatures) Top10"
tissue_type <- "breast"
label <- "human_signatures"
TOP_WHAT <- 10

# Read in ssGSEA and IC50 
load(sprintf("%s/ssGSEA_values_%s_cell_line_ssGSEA_%s_unscaled.RData", output_folder, tissue_type, label))  # created in script 2a
load("IC50_breast.cell.lines_HK.drugs.RData") # created in script 3

# Prepare datasets:
# a) Prepare ssGSEA matrix - cleaning
ssgseaValues <- ssgseaValues[-2,]
ssgseaValues <- as.data.frame(t(ssgseaValues)) #cols = sigs, rows = cell lines
signatures <- colnames(ssgseaValues)
ssgseaValues$Cell_Line <- rownames(ssgseaValues)
ssgseaValues_all <- ssgseaValues

# To remove noise, if "TOP_WHAT" was defined, then we only take the TOP_WHAT (top 5, top 10) cell lines each for BC and LP into account
top_cell_lines <- as.character(top.BC.LP_cell.lines$Cosmic.Identifier)
top_cell_lines <- as.numeric(top_cell_lines)

if(TOP_WHAT != 0){
  # Which cell lines are in the top what
  cell_lines_to_keep <- c(head(top_cell_lines , TOP_WHAT), tail(top_cell_lines , TOP_WHAT))
  # Subset to only these cell lines :D 
  ssgseaValues <- ssgseaValues[which(rownames(ssgseaValues) %in% cell_lines_to_keep),]
}

# b) Prepare IC50 matrix
# Which cell lines  are in cell lines of interest?
rows_to_keep <- which(rownames(IC50_subset) %in% rownames(ssgseaValues))
length(rows_to_keep)
IC50_subset$Cell_Line <- NULL
IC50_subset <- IC50_subset[rows_to_keep,] # rows = cell lines, cols = drugs
drug_identifiers <- colnames(IC50_subset)
IC50_subset$Cell_Line <- rownames(IC50_subset)

# c) Make merged data table of cell lines 
dataset <- merge(x = ssgseaValues, y = IC50_subset, by = "Cell_Line", all.y = F)
IC50_subset$Cell_Line <- NULL
ssgseaValues$Cell_Line <- NULL

# Instantiate empty table
drugCorrelations_all <- data.frame(matrix(ncol=4))
colnames(drugCorrelations_all) <- c("signature", "drug", "pcc", "pvalue")

# For each drug, calculate PCC
for (drugID in drug_identifiers){
  # Perform correlation test
  drugCorrelations = foreach(signature = signatures, .combine = 'rbind') %do% {
    corI = cor.test(dataset[,signature], dataset[,drugID])
    return( data.frame(signature = signature, drug = drugID, pcc = corI$estimate[['cor']], pvalue = corI$p.value))
  }
  # Append as row to main table
  drugCorrelations_all <- rbind(drugCorrelations_all, drugCorrelations)
}
# Remove some R objects
rm(list=c("drugCorrelations", "signature", "drugID"))
# Remove first row
drugCorrelations_all <- drugCorrelations_all[-1, ] 

# save(cell_line_ann, drug_ann, IC50_mat, file="IC50.data.and.annotations.RData")
# load("IC50.data.and.annotations.RData")

# Read in drug annotations
drug_ann <- read.xlsx("TableS1F Drug Ann clean.xlsx")
drug_ann$Name <- make.unique(drug_ann$Name, sep =".")

# Merge drug annotations with correlations
y <- merge(drug_ann, drugCorrelations_all, by.x = "Identifier", by.y = "drug", all.x = F)

# Remove unneccesary columns
colnames(y)[3:8]
drugCorrelations_final <- y[,-c(3:8)]
colnames(drugCorrelations_final)
rm(y)

# library(tidyr)
# pcc_matrix<-spread(drugCorrelations_final[,2:4], key="signature",value="pcc", sep="Name",drop=T)

# Make drug Correlations into a new table, in which each column is a signature
pcc_matrix <- dcast(data = drugCorrelations_final,formula = Name~signature,fun.aggregate = sum,value.var = "pcc")
rownames(pcc_matrix) <- pcc_matrix$Name
pcc_matrix$Name <- NULL

# P Value
pval_matrix <- dcast(data = drugCorrelations_final,formula = Name~signature,fun.aggregate = sum,value.var = "pvalue")
rownames(pval_matrix) <- pcc_matrix$Name
pval_matrix$Name <- NULL

# Create heatmap
pheatmap(pcc_matrix, scale="none", color = colorRampPalette(rev(brewer.pal(8, "RdBu")))(50),
         filename = sprintf("%s/heatmap_GDSC.IC50.PCC_ssgsea_cell.lines_top%s.x.%s_%sfold_unscaled.pdf", output_folder,  TOP_WHAT, label, FC))

# rownames(pcc_matrix_1) <- gsub("\\.1", " (rescreen)", rownames(pcc_matrix_1))

# Rescale data (-1 to 1) and receate heatmap
x_rescale <- rescale(as.matrix(pcc_matrix), c(-1,1))
pheatmap(x_rescale, scale="none", color = colorRampPalette(rev(brewer.pal(8, "RdBu")))(50),
         filename = sprintf("%s/heatmap_GDSC.IC50.PCC_ssgsea_cell.lines_top%s.x.%s_%sfold.pdf", output_folder, TOP_WHAT, label, FC))

# Save R workspace
save.image(sprintf("%s/4a.workspace.RData",output_folder))
