# Addition to Part 2 for HK
output_folder <- "Human Breast Cell Lines IC50 (FC5, p0.05, human_signatures) Top10"
load(sprintf("%s/2.ssGSEA.RData", output_folder))

# Libraries
library(openxlsx)
library(kazutils) # devtools::install_github("kazeera/kazutils")
library(pheatmap)
library(RColorBrewer)

# First rename ssGSEA columns so that the cell line names are shown on the heat map rather than numeric identifiers ----
# Match cell lines from annotation
matches <- match(cell_lines_of_interest, x = colnames(ssGseaValues))
# Check that the matching works
# y <- data.frame(Col = colnames(ssGseaValues), Ann = cell_lines_of_interest[matches], CL =cell_lines_of_interest_name[matches])
# all(y$Col == y$Ann) #T

# Rename columns of ssGSEA to replace COSMIC identifier with Cell Line name
ssGseaValues_named <- ssGseaValues
cell_lines_of_interest_name <- cell_line_ann$Sample.Name[cl_breast_indices]
colnames(ssGseaValues_named) <- cell_lines_of_interest_name[matches]


# Make heat map of enrichment scores
hm_plot <- pheatmap(t(ssGseaValues_named[-2,]),
                 color =  colorRampPalette(rev(brewer.pal(8, "RdBu")))(50), 
                 scale = "column",
                 cellwidth = 50,
                 cellheight = 10,
                 # color = brewer.pal(9, "RdYlGn"),
                 filename = sprintf("%s/heatmap_%s_cell_line_ssGSEA_w.names_%s%s_scaled.by.col_original_premeno.pdf",output_folder, tissue_type, label, FC))

# Rotate dendogram so cluster with higher expressed values are on the top
library(dendextend)
plot(hm_plot$tree_row)

# Rotate
hm_rotated_row_dend <- hm_plot$tree_row %>% rotate(2:1)
# # uncomment the next two lines for what i did for the mouse 3FC top10 signatures
# hm_rotated_row_dend <- hm_rotated_row_dend%>% rotate(4:2)
# hm_rotated_row_dend <- hm_rotated_row_dend%>% rotate(2:1)
# hm_rotated_row_dend <- hm_rotated_row_dend2
# 
# hm_rotated_row_dend <- hm_rotated_row_dend%>% rotate(8:9)
# hm_rotated_row_dend <- hm_rotated_row_dend%>% rotate(2:1)

# Create final heat map
hm_plot_final <- pheatmap(t(ssGseaValues_named[-2,]),
                    color =  colorRampPalette(rev(brewer.pal(8, "RdBu")))(50), 
                    cluster_rows = hm_rotated_row_dend,
                    scale = "column",
                    cellwidth = 50,
                    cellheight = 10,border_color = NA,
                    # color = brewer.pal(9, "RdYlGn"),
                    filename = sprintf("%s/heatmap_%s_cell_line_ssGSEA_w.names_%s%s_scaled.by.col_final_premeno.pdf", output_folder, tissue_type, label, FC))


# Now get differences of top 5 cell lines -------------
df <- data.frame(Cell.Line = colnames(ssGseaValues_named),
                 Cosmic.Identifier = colnames(ssGseaValues),
                 BC.ssGSEA.scores = ssGseaValues_named["BC", ],
                 LP.ssGSEA.scores = ssGseaValues_named["LP", ])
df$Difference.BC.minus.LP <- df$BC.ssGSEA.scores - df$LP.ssGSEA.scores
# Sort highest to lowest difference
df <- df[order(-df$Difference.BC.minus.LP),]#- means we're sorting in descending order

# Label Top 5 BC and Top 5 LP cell lines
tail(df$Difference.BC.minus.LP, n = 5)
df$Top.5 <- c(rep("Top.5.BC", 5), rep("", times = nrow(df)-10), rep("Top.5.LP", 5))
df$Top.10 <- c(rep("Top.10.BC", 10), rep("", times = nrow(df)-20), rep("Top.10.LP", 10))
top.BC.LP_cell.lines <- df

# Save data to Excel file
wb <-createWorkbook()
write_list_to_Excel(sheet_name = sprintf("Human Signatures FC>%s pval<%s", FC, pval_thres),
                    my_list = signatures, wb = wb)
write_df_to_Excel(sheet_name = "ssGSEA Values", my_df = t(ssGseaValues_named), wb, incl_rownames = TRUE)
write_df_to_Excel(sheet_name = "Top 5 Scores in BC and LP", my_df = df, wb, incl_rownames = FALSE)
excel_filename <- sprintf("%s/%s_ssGSEA_premenoHumanMammSigFC%s_x_GDSCcellLine.xlsx", output_folder, format(Sys.Date(), "%y%m%d"), FC)
saveWorkbook(wb, file = excel_filename,overwrite = T) # excel filename can't be too long

# Save image
save.image(sprintf("%s/2b.workspace.RData", output_folder))
  