# Run this script to get a heat map from the protein expression matrix
# All expression values are LFQ-adjusted IBAQ values (proteome).	
# The 0s in the data were imputed with a random number between 1 and 1.5.	

#helper functions
get_present_gene_name <- function(gene_column, signatures){
  # returns a list of genes, where the position of the gene name matches the position in the original gene_list (in which gene names are concatenated by ";")
  # if there is a gene name, it is present in the signatures, if not then it is absent
  # e.g. INPUTS gene_list = "a;b;c", "c;d;f", "d", "l;d;a"  
  #             signatures = "l", "q", "b     
  #      OUTPUT "BC", "", "", "l" 
  
  # Split the gene symbols delimited by ";" into seperate entities as elements of a list
  gene_column <- strsplit(gene_column,";")
  # Remove the ".number" that was appended to make the gene names unique
  gene_column <- lapply(gene_column, function(x) gsub("\\..*","",x))
  # Look for any genes present in the signature in each split list vector 
  present_genes <- lapply(gene_column, function(x){
    if(any(x %in% signatures)){
      genes <- x[x %in% signatures]  #find whether genes are in signatures
      genes <- paste(genes, collapse=";") #concatenate using ";" if necessary
    } 
    else return("")
  })
  # unlist(sapply(1:6222, function(x){ if(length(present_genes[[x]]) > 1) return(paste(present_genes[[x]], collapse=";")) })
  # Unlist list to vector form
  present_genes <- unlist(present_genes)
  
  # Return a vector of genes that are present in the signatures ("" if absent), conserving the order of the original gene list
  return(present_genes)
}

subset_proteome_table2 <- function(table, gene_column, signatures){
  #returns a subset of a proteome "table" with only the signatures 
  #subset original proteome data based on genes that match the genes of interest list
  present_genes <- get_present_gene_name(gene_column, signatures)
  table <- table[present_genes != "",]
  rownames(table) <- make.unique(present_genes[present_genes != ""], sep = ".")
  return (table)
}


# Arguments: 
# full proteome - all protein expression data, matrix where row names correspond to proteins/gene symbols and columns correspond to sample info
# signatures - list of genes to subset
# title - name of heat map and filename id
# celltypes - vector with (B, ML, LP, S cell types) according to order of full proteome
# average - logical. If true, will map the average expression values across replicates for each cell type
make_human_heatmap <- function(full_proteome, signatures, title, average){
  # Call libaries (usually at beginning of script)
  library(pheatmap)
  library(RColorBrewer)
  library(cluster)
  library(scales)
  # Define variables for heat map ----
  # Make annotation column (sample information)
  # celltypes <- c("BC", "ML", "LP")
  # num_replicates <- 10  #9 for basal cell
  # phases <- c("Follicular", "Luteal", "PostMenopausal")
  
  # Subset proteome ----
  raw_value_mat <- subset_proteome_table2(table = full_proteome, gene_column = rownames(full_proteome), signatures)
  
  # Compute z-scores to plot on heat map
  z_score_mat <- apply(raw_value_mat, 1, function(x) (x - mean(x)) / sd(x))
  exp_matrix <- t(z_score_mat)
  exp_matrix <- rescale(as.matrix(exp_matrix), to = c(-2, 2))
  
  # Make hierarchical cluster objects (dendogram - like structures)----
  # Clustering method ----
  cluster_method <- "pearson" 
  
  # b) Cluster rows (genes)
  hc_proteins <- as.dist(1-cor(t(exp_matrix), method =cluster_method))
  #use = "pairwise.complete.obs",
  hc_proteins <- as.dendrogram(diana(hc_proteins)) #DIvisive ANAlysis Clustering
  
  
  # Make annotations -----------------    
  #make annotation column (sample information)
  annotation_col <- NA
  annotation_col <- data.frame(
    CellType = substr(colnames(exp_matrix_z), 14, 15)
    # CellType = factor(rep(celltypes, each = num_replicates))
    # PatientAge = rep(c("31-40","31-40","31-40", "<30","<30","31-40", "51-60", "51-60",  ">60", ">60"), length(celltypes)),
    # Phase = rep(rep(phases, c(3,3,4)), length(celltypes))
  )
  # # Drop 9th row (patient sample BC_49_15)
  # annotation_col <- annotation_col[-9,] 
  # annotation_col <- data.frame(annotation_col)
  rownames(annotation_col) <- colnames(exp_matrix)
  # 
  # Make annotation row (gene information) -- not used right now
  annotation_row = NA
  # all_genes <- rownames(exp_matrix)
  # rownames(annotation_row) <- rownames(exp_matrix)
  
  # # Define annotation colors 
  # # a) Pick colors for the age groups 
  # newAgeCols <- colorRampPalette(c("white", "black")) #"yellowgreen
  # myAgeColors <- newAgeCols(length(unique(annotation_col$PatientAge))+1)
  # names(myAgeColors) <-c("<30", "31-40", "41-50", "51-60",  ">60")
  # 
  # # b) Pick colors for the age groups 
  # newPhaseCols <- colorRampPalette(c("lightgoldenrodyellow", "goldenrod2"))
  # myPhaseColors <- newPhaseCols(length(phases))
  # names(myPhaseColors) <- phases
  # 
  # c) Pick colors for the cell-types
  myCellColors = c(BC="red", LP="blue", LC="darkgreen")#, S="darkgray")
  
  # d) Make a list of the annotation colors
  annotation_colors <- list(
    CellType = myCellColors
    # Phase = myPhaseColors,
    # PatientAge = myAgeColors 
  )
  
  # Make heat map ----------
  # a) Define title and filename if desired
  title_for_plot <- sprintf("%s (%s proteins)", title, nrow(exp_matrix))
  filename <- paste(format(Sys.Date(), "%y%m%d"), title, "human_mamm_heatmap.jpg", sep="_")# "human_proteome2.pdf" #pheatmap can save to multiple file types
  filename <- sprintf("%s/%s", out_dir, filename)
  # b) pheatmap function call
  pheatmap(exp_matrix, 
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           #color = colorRampPalette(c("chartreuse2", "white", "darkorchid4"))(50),
           color = colorRampPalette(c("darkolivegreen", "darkgreen", "white", "darkorchid3", "darkorchid4"))(50), #deeppink4
           # color = bluered(75),
           scale = "none", #row",
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           cluster_cols = F,#hc_samples,
           cluster_rows = as.hclust(hc_proteins),
           # clustering_method = hc_proteins$dist.method,
           border_color = NA,
           cellheight = 3,
           # cutree_rows = 11,
           cellwidth = 15,
           show_rownames = TRUE,
           show_colnames = TRUE,
           treeheight_col = 10,
           treeheight_row = 50,
           fontsize_row = 3,
           main=title_for_plot,
           filename = filename
  )
  #-----------------   
}
