# Colors to annotate BC subtypes
PAM50_colors <- c(LumA="mediumpurple1", LumB="mediumpurple4", Basal="seagreen1", Normal="seagreen", HER2="goldenrod1")


# Heatmap function
make_heatmap <- function(mat, cols = NA, label2="", scale="column", pal="RdBu", BC_subtypes = F, ann_row = NA, ann_col = NA, ann_colors = F){
  library(RColorBrewer)
  library(pheatmap)
  
  
  if(BC_subtypes){
    # Get annotations, for samples (rows)
    # col_ann$final_processedID == colnames(ssgsea_val)
    ann_row <- data.frame(col_ann[,c(cols)], 
                          row.names = col_ann$final_processedID)
    colnames(ann_row) <- cols
  }
  
  # Annotation colors
  cols <- colnames(ann_row)
  ann_colors <- lapply(cols, function(x){
    y <- ann_row[,x] %>% as.character %>% unique
    
    if(any(y%in%names(PAM50_colors))){
      PAM50_colors[names(PAM50_colors) %in% y]  
      
    }else{
      get_col_palette(ann_row[,x], brew_pal = "Accent")
    }
  })
  names(ann_colors) <- cols
  
  # Make heat map of enrichment scores
  pheatmap::pheatmap(mat,
           color = colorRampPalette(rev(brewer.pal(8, pal)))(50), 
           annotation_row = ann_row,
           annotation_colors = ann_colors,
           annotation_col = ann_col,
           scale = scale,
           cellwidth = 50,
           cellheight = 10,
           main = label2,
           # color = brewer.pal(9, "RdYlGn"),
           filename = sprintf("%s/heatmap_scaled_ann_%s.pdf",output_folder, label2))
  
}

# Get colours for a vector #got this from Mathepan's heatmap function
# v = unique elements
#' custom - vector with 2 values low_col and high_col are colors in RColorBrewer that are two ends of a gradient 
#' color_pallette - specific pallette from RColorBrewer
get_col_palette <- function(v, brew_pal, custom=NA, rev=F, rearr=F){
  require(RColorBrewer)
  # Get unique elements and rearrange
  v <- unique(v)
  # # Rearrange   # 
  if(rearr)
    v <- sample(v)
  
  if(!is.na(custom)) {
    # Make color pallette with number of elements in v
    newCols <- colorRampPalette(custom)#c(low_col, high_col))
    myColors <- newCols(length(v))
  }
  if(!is.na(brew_pal)){
    # Get max colors of palette
    max_n <- brewer.pal.info$maxcolors[grep(brew_pal, rownames(brewer.pal.info))]
    # Create brewer pal
    # Reverse if required
    if(rev){
      brew <- rev(brewer.pal(n = max_n, name = brew_pal))
    }else{
      brew <- brewer.pal(n = max_n, name = brew_pal)
    }
    # Finally, create gradient
    myColors <- colorRampPalette(brew)(length(v))
  }
  names(myColors) <- v
  return(myColors)
}