#' Purpose:
#' To find finds enrichment of damage genes in all cells and GMNN+ cells.
#' 
#' Added cell cycle phase annotation in script 1.R
#' Added "RAD51_set", "HR", "NHEJ" in script 2.R

# Remove all objects in environment except b2 (b2 = Seurat object)
rm(list=ls()[!ls()%in%"b2"])

# Import R libraries
library(kazutils)
load_packages(c("Seurat", "ggplot2", "cowplot", "dplyr"))

# Import functions
source("../plot_boxplot.R")

# Out directory
out_dir <- kazutils::create_folder("Plots")

# Data 
b2 <- readRDS("../annotated_Seurat_b2.rds")

# Colors
cell_colors <- readRDS("../colors_broadcelltype2.rds")
colors_Phase <- readRDS("../colors_Revelio_ccPhase.rds")

# Visualize enrichment scores
sets <- c("RAD51_set", "HR", "NHEJ")
ES2 <- data.frame(b2[[]], Idents(b2))
ES2$broadcelltype2 <- factor(ES2$broadcelltype2, levels = c("LP", "ML", "BC"))

# Only look at HR and RAD51 in GMNN positive cells for HR and RAD51
b2_GMNN <- subset(b2, GMNN > 0)
ES3 <- data.frame(b2_GMNN[[]], Idents(b2_GMNN))

cell_number <- table(ES3$broadcelltype2) %>% paste(names(.), ., collapse=".")

for (set in c("RAD51_set", "HR")){
  # Data frame
  df <- ES3[,c("Idents.b2_GMNN.",set)]
  # Make plot
  p1 <- plot_boxplot(df, colors = cell_colors, mean_line=T, ylab="Enrichment Score", 
                     title = sprintf("Enrichment of %s", set), subtitle = paste(cell_number, "in GMNN+ cells. Line=mean", sep=";"))
  p2 <- plot_boxplot(df, colors = cell_colors, box = T, ylab="Enrichment Score",
                     title = sprintf("Enrichment of %s", set), subtitle = paste(cell_number, "in GMNN+ cells", sep=";"), add_Tukey = F)
  
  # Plot
  pdf(sprintf("%s/3_onlyGMNNpositive_%s.pdf", out_dir, set))
  print(p1)
  print(p2)
  dev.off()
  
  # Make boxplot for each cell type
  plots <- lapply(unique(b2$broadcelltype2), function(celltype){
    # Data frame
    df <- ES3[ES3$broadcelltype2 == celltype, c("Phase",set)]
    
    # Make plot
    plot_boxplot(df, colors = colors_Phase, mean_line=T, box = T,
                 title = sprintf("Enrichment of %s in %s", set, celltype), subtitle = "GMNN+ cells")
  })
}

