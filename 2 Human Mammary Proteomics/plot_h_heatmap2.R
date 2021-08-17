# plot (human) heatmap function takes in the following to create a heat map..
# Arguments:
# exp_matrix = numeric matrix of expression data with columns as samples and genes as rows
# hc_samples, hc_proteins = hc_clust objects indicating hierarchal clustering 
# num_replicates = number of replicates of each cell type (synonymous with number of patients)
# genes_of_interest = vector of genes used to annotate rows
# genes_of_interest_label = string indicating group name of these genes of interest
# filename =  string - name of the file with directory and format as .pdf, .jpeg, etc
# title = string - title of the heat map 
# celltypes = vector of the cell type names eg. c("B", "LM", "LP", "S")

# heat map colors:
# purple - low 
# green - high 

plot_h_heatmap_premeno <-function(exp_matrix, hc_samples, hc_proteins, filename, title, pint.pheno, cluster_vector){
  # Call libaries (usually at beginning of script)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)

  exp_matrix <- rescale(exp_matrix, to = c(-2, 2))

  
  # celltypes <- c("BC", "LM", "LP")
  # num_replicates <- 10  #9 for basal cell
  phases <- c("Follicular", "Luteal")
  
  # Make annotations -----------------    
  #make annotation column (sample information)
  annotation_col <- data.frame(
    CellType = gsub("LC", "LM", pint.pheno$cell_type),
    # PatientAge = rep(c("31-40","31-40","31-40", "<30","<30","31-40", "51-60", "51-60",  ">60", ">60"), length(celltypes)),
    PatientAge = pint.pheno$patient_age_group,
    HormoneStatus = pint.pheno$hormone_status
  )
  rownames(annotation_col) <- colnames(exp_matrix)
  
  
  #make annotation row (gene information)
  annotation_row = data.frame(
    SignificantHits = cluster_vector
  )
  rownames(annotation_row) <- rownames(exp_matrix)
  
  # Define annotation colors -----------------
  newAgeCols <- colorRampPalette(c("white", "black")) #"yellowgreen
  myAgeColors <- newAgeCols(length(unique(annotation_col$PatientAge))+1)
  names(myAgeColors) <- unique(annotation_col$PatientAge)
  
  newPhaseCols <- colorRampPalette(c("lightgoldenrodyellow", "goldenrod2"))
  myPhaseColors <- newPhaseCols(length(phases))
  names(myPhaseColors) <- phases
  
  annotation_colors <- list(
    CellType = c(BC="red", LP="darkgreen", LM="blue"),#, S="darkgray"),
    HormoneStatus = myPhaseColors,
    PatientAge = myAgeColors,
    SignificantHits =  c(BC="red", LP="blue", LM="darkgreen", Unassigned="white")
    #Pathway = c(genes[genes%in%hdr_genes]="blue", genes[!genes%in%hdr_genes]="black")
  )
  # filename = sprintf("%s_possemato_heat_map.png", format(Sys.Date(), "%Y%m%d"))
  # title = NA
  
  # pheatmap function call -----------------    
  pheatmap(exp_matrix, 
           #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           #color = colorRampPalette(c("chartreuse2", "white", "darkorchid4"))(50),
           # color = colorRampPalette(c("darkorchid3", "darkorchid4", "white", "forestgreen", "darkolivegreen"))(50), #deeppink4
           color = colorRampPalette(c("darkgreen", "forestgreen", "white", "darkorchid3", "darkorchid4"))(50),
           # #color = bluered(75),
           scale = "none", #"row",
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           # cluster_cols = hc_samples,
           # cluster_rows = hc_proteins,
           # clustering_method = hc_proteins$dist.method,
           border_color = NA,
           # cutree_rows = 10,
           cellwidth = 15,
           show_rownames = FALSE,
           show_colnames = FALSE,
           treeheight_col = 10,
           treeheight_row = 50,
           # cellheight = 1,
           # fontsize_row = 1,
           # height = 10,
           main=title,
           filename = filename
  )
  #-----------------    
  # accepted by pheatmap: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
}
