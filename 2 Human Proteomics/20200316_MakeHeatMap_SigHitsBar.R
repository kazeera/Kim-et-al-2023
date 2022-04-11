#' Purpose:
#' To create heatmaps for each signature and a final heat map with "cluster assignment" bar
#' Input: 
#'  gene lists
#' Output:
#'  heat maps - individual for each signature, final heat map with "cluster assignment" bar
#'  RData file
#'  Final heat map
#' 
#' 20200316: Edited script 4c from https://github.com/mcclo/Mahendralingam-et-al.-Nat-Metab/ for Hyeyeon's premenopausal total proteome human signature
#' Edits:
#' - made input files more general (lapply, RData)
#' - replaced "ML" occurences with "LM" - in plot_heatmap2.R too
#' - made a function to make dendograms
#' 
# Note: significant_hits object holds the final signatures
# BC = basal cell, LM = luminal mature, LP = luminal progenitor

# Clear R environment
rm(list=ls())

# Call required packages into environment
library(scales)
library(pheatmap)
library(RColorBrewer)
library(cluster) # for make_human_heatmap.R

# Load files from previous scripts
lapply(list.files(pattern = "RData"), load, .GlobalEnv)
rm(signatures)

# Create output folder
out_dir <- sprintf("%s_premeno_hm_%sFC", format(Sys.Date(), "%Y%m%d"), FC)
dir.create(out_dir)

# computes and returns matrix with z scores for heatmap plotting
compute_z_scores <- function(matrix){
  # Formula for z score: 
  #         [element x - mean of row x is in]/ standard deviation of row x is in
  return (apply(matrix, 1, function(x) (x - mean(x)) / sd(x)))
}
# Make z score matrix
exp_matrix <- h.pint.combat_matrix_no.post.meno
rm(h.pint.combat_matrix_no.post.meno)
exp_matrix_z <- t(compute_z_scores(exp_matrix))

# # Make dendograms
# # load package for diana function
# library(cluster)
# library(gplots)
# sample_dist <- as.dist(1-cor(exp_matrix, method = "euclidean"))
# sample_dend <- as.dendrogram(diana(sample_dist))
# 
# protein_dist <- as.dist(1-cor(t(exp_matrix), method = "pearson"))
# # prot_dend <- as.dendrogram(diana(sample.dend)) 

# Make a new vector that indicates which genes belong to which cell type
significant_hits_bar <- ifelse(rownames(exp_matrix) %in% significant_hits$BC, "BC",
                                      ifelse(rownames(exp_matrix) %in% significant_hits$LP, "LP",
                                             ifelse(rownames(exp_matrix) %in% significant_hits$LM, "LM", "Unassigned")))

# Define order
celltype_gene_assignments <- c("BC", "LM", "LP", "Unassigned")

# Plot final heatmap with cluster vector
source("plot_h_heatmap2.R")
filename <- sprintf("%s_FINAL_premeno_total_heat_map_signatures.%sFC.pdf", format(Sys.Date(), "%Y%m%d"), FC)
filename <- sprintf("%s/%s", out_dir, filename)

final_heatmap_plot <- plot_h_heatmap_premeno(exp_matrix = exp_matrix_z, hc_samples = T, 
                                      hc_proteins = T, filename = filename, title = NA,
                                      pint.pheno = pint.pheno_no.meno, cluster_vector = significant_hits_bar)

# Make heatmap for each cell type
source("make_human_heatmap.R")
# Define function to make heat maps for each of the 3 cell signatures
#' @param cell_signatures is a list of 3 signatures for basal, LM and LP respectively
make_cell_signature_heatmap <- function(cell_signatures, title){
  # Names of titles/file names
  plot_titles <- c("BC_population", "LM_population", "LP_population") 
  celltypes <- c("BC", "LM", "LP")
  # Iteratively make heatmap
  sapply(1:3, function(i){
    # Make title
    title <- paste(plot_titles[[i]], title, sep="_")
    # Get gene list
    signature <- cell_signatures[[i]]
    # Subset and create heatmap
    x <- make_human_heatmap(exp_matrix, signature, title, average = FALSE)
    # dev.off()
  })
}

# Run function for each signature list (before, after ANOVA, and after ANOVA+Tukey)
make_cell_signature_heatmap(significant_hits, "significant_hits") #unclustered #cluster_cols = F






