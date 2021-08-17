# creates heat maps by calling helper functions 
plot_genes_of_interest <- function(genes_of_interest, name){
  # # Make an output directory
  # output_dir <- sprintf("%s mouse heat map", name)
  # dir.create(output_dir)
  # 
  # CREATE Z SCORE MATRIX -----------------
  # Subset proteome z scores with genes of interest
  m_exp_matrix <- subset_proteome_table3(table = m.z.scores, 
                                         gene_column = m.z.scores$PROBE.ID, 
                                         signatures = genes_of_interest)
  colnames(m_exp_matrix)
  m_exp_matrix_original <- m_exp_matrix
  # Drop probe ID column
  m_exp_matrix$PROBE.ID <- NULL
  # m_exp_matrix -> m_exp_matrix_all
  m_exp_matrix <- m_exp_matrix[,grep("EP",colnames(m_exp_matrix))]
  colnames(m_exp_matrix)
  
  # Define columns annotations (sample descriptions)
  cell_type_samples <- rep(c("B", "LP", "LM"), each = 2)
  hormone_samples <- rep("EP", times = 6)
  replicate <- rep(1:2, times = 3)
  
  
  # --------------------
  
  # Make dendograms
  dist_samples <- as.dist(1-cor(m_exp_matrix,
                                #use = "pairwise.complete.obs",
                                method = "spearman"))
  hclust_samples <- hclust(dist_samples)
  # png(sprintf("sample_dendogram_%s.png", name))
  # plot(as.dendrogram(hclust_samples)) 
  # dev.off()
  ## NB: can combine into one line of code but gets confusing
  
  # Protein dendrogram m.pint
  dist_proteins <- as.dist(1-cor(t(m_exp_matrix),
                                 #use = "pairwise.complete.obs",
                                 method = "spearman"))
  # be prepared - this next line of code takes like 5 minutes
  hclust_proteins <- hclust(dist_proteins) 
  ###
  # save(m_exp_matrix, hclust_samples, dist_samples, hclust_proteins, dist_proteins, file = "mouse_DDR.pw_dendograms.RData")
  #-------------------------------------------------------
  title_for_plot <- sprintf("%s of %s genes found",  nrow(m_exp_matrix), length(genes_of_interest))
  filename <- sprintf("%s_HK_mouse.proteomics_DDR.pw.from.cytoscape_dist.spearman_linkage.complete_heatmap_%s.pdf", format(Sys.Date(), "%y%m%d"), name)
  
  plot_mouse_heatmap(m_exp_matrix, filename, cell_type_samples, hormone_samples,hclust_proteins, hclust_samples, title_for_plot)
  return(m_exp_matrix_original)
}


#returns a subset of a proteome "table" with only the signatures 
#subset original proteome data based on genes that match the genes of interest list
subset_proteome_table3 <- function(table, gene_column, signatures){
  subsett<- gene_column%in%signatures
  table2 <- table[subsett,]
  rownames(table2) <- rownames(table)[subsett]
  return (table2)
}


plot_mouse_heatmap <- function(exp_matrix, filename, cell_type_samples, hormone_samples, hc_proteins, hc_samples, title_for_plot){
  library(pheatmap)
  library(RColorBrewer)
  # Make hierarchical cluster objects (dendogram - like structures)
  # a) Define clustering method 
  # possible methods: "euclidean", "maximum", "manhattan", "canberra","minkowski"
  # cluster_method <- "euclidean" #this is the default
  # # b) Cluster rows (samples)
  # hc_samples <- hclust(dist(exp_matrix, method=cluster_method))
  # 
  # # c) Transpose the exp_matrix and cluster columns(genes)
  # hc_proteins <- hclust(dist(t(exp_matrix), method=cluster_method))
  # 
  rownames(exp_matrix) <- gsub(";.*", "", rownames(exp_matrix))
  # Make annotation column (sample information)
  annotation_col <- data.frame(
    CellType = cell_type_samples
    # Hormone = hormone_samples
  )
  # annotation_row <- data.frame(
  #   CellType = rep(c("B",  "LP", "LM"), each = 4),
  #   Hormone = rep(c("E", "EP"), each = 2, times = 3)
  # )
  rownames(annotation_col) <- colnames(exp_matrix)
  
  # Make annotation column (gene information) -- not used right now
  annotation_row <- NA 
  
  # Make a list of the annotation colors
  annotation_colors <- list(
    CellType = c(B="red", LP="skyblue", LM="darkgreen")
    # Hormone = c(E="bisque", EP="black")
  )
  
  # Make heat map 
  # a) Define title and filename if desired
  # filename <- paste(format(Sys.Date(), "%y%m%d"), title, "mouse_heatmap.pdf", sep="_")# "human_proteome2.pdf" #pheatmap can save to multiple file types
  
  # b) pheatmap function call
  pheatmap(exp_matrix, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           # color = colorRampPalette(c("chartreuse2", "white", "darkorchid4"))(50),
           # color = colorRampPalette(c("darkolivegreen", "darkgreen", "white", "darkorchid3", "darkorchid4"))(50), #deeppink4
           # color = colorRampPalette(c("darkgreen", "white", "darkorchid3"))(100), #deeppink4
           # color = bluered(75),
           scale = "row",
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_cols = hc_samples,
           cluster_rows = hc_proteins,
           clustering_method = "complete",
           border_color = NA,
           # cellheight = 10,
           cellheight = 6,
           # cutree_rows = 11,  
           cellwidth = 10,
           show_rownames = TRUE,
           show_colnames = F,#TRUE,
           treeheight_col = 15,
           treeheight_row = 50,
           fontsize_row = 6,
           fontsize_col = 1,
           main=title_for_plot,
           filename = filename
  )
  
}


# This function writes a data frame with column names to Excel
# It will make a worksheet in the workbook (wb) specified with the name of "current_sheet"
write_df_to_Excel <- function(current_sheet, my_df, wb, incl_rowNames){
  # add a worksheet to the workbook with the name specified in the argument 
  addWorksheet(wb, current_sheet)
  # Write data frame to sheet
  writeData(wb, x = my_df, sheet = current_sheet, colNames = T, rowNames = incl_rowNames)
  return(wb)
}
