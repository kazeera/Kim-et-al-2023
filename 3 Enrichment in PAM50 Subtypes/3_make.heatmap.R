# Run this script to get a heat map from the ssGSEA subtype matrix
# rm(list = ls())
#Prepare matrix
exp_matrix <- t(ssGSEA.subtype[,2:4]) 
colnames(exp_matrix) <- rownames(ssGSEA.subtype[,2:4])
rownames(exp_matrix) <- colnames(ssGSEA.subtype[,2:4])
sample_subtypes <- ssGSEA.subtype$subtype

# Call libaries (usually at beginning of script)
library(pheatmap)
library(RColorBrewer)
library(ggsci)

# Arguments: 
# exp_matrix -rows = name of cell type, columns = ssGSEA correlation
# signatures - list of genes to subset
# title - name of heat map and filename id
# subpopulations - vector with (B, ML, LP, S cell types) according to order of full proteome
make_ssGSEA_heatmap <- function(){
  sample_subtypes <- ssGSEA.subtype$subtype
  subtypes <- c("Basal", "claudin-low", "Her2", "LumA", "LumB") 
  title = "Celltypes_cancer_subtypes"
  subpopulations <- colnames(exp_matrix)
  
  # Make hierarchical cluster objects (dendogram - like structures)
  # a) Define clustering method 
  # possible methods: "euclidean", "maximum", "manhattan", "canberra","minkowski"
  cluster_method <- "euclidean" #this is the default
  linkage_method <- "complete"
  # b) Cluster rows (genes)
  hc_rows <- hclust(dist(exp_matrix, method=cluster_method), method = linkage_method)
  # c) Transpose the exp_matrix and cluster columns (samples)
  hc_cols <- hclust(dist(t(exp_matrix), method=cluster_method), method = linkage_method)
  
  # Make annotation column (sample information)
  annotation_col <- data.frame(Subtype = sample_subtypes)
  rownames(annotation_col) <- colnames(exp_matrix)
  
  # Make annotation row (samples)
  annotation_row <- data.frame(CellType = c("B", "ML", "LP"))
  rownames(annotation_row) <- rownames(exp_matrix)
  
  # Define annotation colors 
  # # a) Pick colors for the age groups 
  # mySubtypeColours <- brewer.pal(n = 5, name = "Dark2")
  # names(mySubtypeColours) <- subtypes
  
  mySubtypeColours <- ggsci::pal_jco("default")(5)
  names(mySubtypeColours) <- subtypes
  # b) Pick colors for the cell-types
  myCellColors = c(B="red", LP="skyblue", ML="blue")
  
  # c) Make a list of the annotation colors
  annotation_colors <- list(
    CellType = myCellColors,
    Subtype = mySubtypeColours
  )
  
  # 14. Make heat map 
  # a) Define title and filename if desired
  filename <- sprintf("%s/%s_ssGSEA_human.mamm.sig_METRABRIC.BC.subtypes_heatmap.pdf",output_folder,todays_date)#" #pheatmap can save to multiple file types
  
  # b) pheatmap function call
  pheatmap(exp_matrix, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           # color = colorRampPalette(c("chartreuse2", "white", "darkorchid4"))(50),
           # # color = bluered(75),
           scale = "column",
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_cols = hc_cols,
           cluster_rows = hc_rows,
           # clustering_method = "ward.D2",
           border_color = NA,
           cellheight = 15,
           # cutree_rows = 11,
           # cellwidth = 1,
           show_rownames = TRUE,
           show_colnames = FALSE,
           # treeheight_col = 10,
           # treeheight_row = 50,
           fontsize = 5,
           main=title_for_plot,
           filename = filename
  )
    #-----------------   
}
title_for_plot <- signature_file
make_ssGSEA_heatmap()
