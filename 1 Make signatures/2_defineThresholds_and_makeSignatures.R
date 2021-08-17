# continuation of script
# call library
library(openxlsx)

celltypes <- c("BC", "LM", "LP")

# Make output folder
output_folder <- sprintf("Human Mammary Signatures (FC%s, p-val%s)", FC, pval_thres)
dir.create(output_folder)

sig_dfs <- lapply(1:3, function(i){
  # add column for significance - value will be "Yes" if pass both thresholds
  add_sig_columns(subset_dfs[[i]], pval_thres, log2FC_thres_up)
})
names(sig_dfs) <- celltypes
names(subset_dfs) <- celltypes

# 3 tabs -  Write important data frames and lists to Excel file
excel_filename <- sprintf("%s/%s_tables_ANOVA_Tukey_p.%s_FCthres.%s_premeno.only.xlsx", output_folder, format(Sys.Date(), "%Y%m%d"), pval_thres, FC)
wb <- openxlsx::createWorkbook("signatures")

#
significant_hits <- vector("list", 3)
names(significant_hits) <- celltypes

for(i in 1:3){
  df <- sig_dfs[[i]]
  celltype <- names(sig_dfs[i])
  
  x <- sprintf("Analyzing %s..", celltype)
  print(x)
  
  current_sheet <- sprintf("%s Significance Table", celltype)
  addWorksheet(wb, current_sheet)
  print(current_sheet)
  
  description1 <- sprintf("conditions: p_val < %s, log2FC >= %s", pval_thres, log2FC_thres_up)
  description2 <- sprintf("%s proteins meet condition for %s", length(which( df$Both_Conditions_Met == "Yes")), celltype)
  hits <- data.frame(SignificantHits = df$Gene[which(df$Both_Conditions_Met == "Yes")])
  
  writeData(wb, sheet = current_sheet, x = description1, startRow = 1)
  writeData(wb, sheet = current_sheet, x = description2, startRow = 2)
  writeData(wb, sheet = current_sheet, x = df, startRow = 4)
  
  writeData(wb, sheet = current_sheet, x = hits, startRow = 4, startCol = ncol(df)+5)
  significant_hits[[i]] <- as.character(hits$SignificantHits)
}
significant_hits <- lapply(significant_hits, FUN = sort)

# log2FC
addWorksheet(wb, "All Tukey and log2FC Results")
writeData(wb, sheet = "All Tukey and log2FC Results", x = cbind(significance_df, log2FC_df), rowNames = FALSE)

# clean signatures for Enrichr
signatures <- lapply(significant_hits, function(x){
  x <- gsub(";.*", "", x) #remove any alias after semicolons
  x <- gsub("\\..*","",x) #remove any period and number appended to make gene name unique
  x <- x[x!=""]
  x <- unique(x)
  return(x)
})
write_list_to_Excel("Signatures (Unique Genes)", signatures, wb)
write_list_to_Excel("Signatures", my_list = significant_hits, wb)

# Combat proteome
addWorksheet(wb, "Total Proteome Combat")
writeData(wb, sheet = "Total Proteome Combat", x = h.pint.combat, rowNames = TRUE)

#Save to excel
saveWorkbook(wb, file = excel_filename, overwrite = TRUE)

# Print to console
# get lengths of signatures
for(i in 1:3){
  x <- sprintf("Length of %s significant hits: %d", 
               names(sig_dfs)[i], length(significant_hits[[i]]))
  print(x)
  x <- sprintf("Length of %s signature (unique gene symbols): %d", 
               names(sig_dfs)[i], length(signatures[[i]]))
  print(x)
}

#--------------------------------------
# File 2 for Signatures files - table of signatures - needed for Luis' ssGSEA analysis
excel_filename2 <- sprintf("%s/3-signatures_FC%s_premeno.only.xlsx", output_folder, FC)
cell_names <- c("Basal", "Luminal Mature", "Luminal Progenitor")
significant_hits_named <- signatures
names(significant_hits_named) <- cell_names

# Create workbook and write lists to it
wb2 <- createWorkbook()
write_list_to_Excel("Signatures", my_list = significant_hits_named, wb2)
print(sprintf("Saving Excel file with signatures only (%s)..", excel_filename2))
saveWorkbook(wb2, file = excel_filename2, overwrite = TRUE)

#-------------- save RData for making heatmap
RData_filename <- sprintf("%s/human.signatures_FC%s.pval%s_premeno.only.RData", output_folder, FC, pval_thres)
print(sprintf("Saving RData file with signatures (%s)..", RData_filename))
save(significant_hits, signatures, FC, log2FC_thres_up, pval_thres, file=RData_filename)
