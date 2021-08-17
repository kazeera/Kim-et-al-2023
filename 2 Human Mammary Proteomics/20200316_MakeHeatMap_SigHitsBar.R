# Part 2 computes the final signatures - after ANOVA and Tukey test (p<0.05)

# Note: significant_hits object holds the final signatures
# BC = basal cell, LM = luminal mature, LP = luminal progenitor

# Input: gene lists
# Output:
#' heat maps - individual for each signature, final heat map with "cluster assignment" bar
#' RData file
#' Final heat map
 
#' 20200316 edit: Edited script from Mathepan's folder: D:\KA\Mathepan\Mathepan - Part 3 With New Matched Data\Analysis included in First Submission to Nat Met\1E Metabolic Heat Map and Signatures\2 Make Heat Map
#' called "20190607_MakeHeatMap_SigHitsBar.R" for Hyeyeon's premenopausal total proteome human signature
#' edits:
#' - made input files more general (lapply, RData)
#' - replaced "ML" occurences with "LM" - in plot_heatmap2.R too
#' - made a function to make dendograms
#' - have an out_directory #TODO figure out options() to redirect output in R

rm(list=ls())
# Call required packages into environment
library(pheatmap)
library(RColorBrewer)

# load files from previous scripts
lapply(list.files(pattern = "RData"), load, .GlobalEnv)
rm(signatures)

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

# define order
celltype_gene_assignments <- c("BC", "LM", "LP", "Unassigned")


# 7. Plot final heatmap with cluster vector----------------------------------
source("plot_h_heatmap2.R")
filename <- sprintf("%s_FINAL_premeno_total_heat_map_signatures.%sFC.pdf", format(Sys.Date(), "%Y%m%d"), FC)
filename <- sprintf("%s/%s", out_dir, filename)

final_heatmap_plot <- plot_h_heatmap_premeno(exp_matrix = exp_matrix_z, hc_samples = T, 
                                      hc_proteins = T, filename = filename, title = NA,
                                      pint.pheno = pint.pheno_no.meno, cluster_vector = significant_hits_bar)

# a) Make heatmap for each 
source("make_human_heatmap.R")
# Define function to make heat maps for each o the 3 cell signatures-----------------
make_cell_signature_heatmap <- function(cell_signatures, title){
  # cell_signatures is a list of 3 signatures for basal, LM and LP respectively
  plot_titles <- c("BC_population", "LM_population", "LP_population") #names of titles/file names
  celltypes <- c("BC", "LM", "LP")
  sapply(1:3, function(i){
    title <- paste(plot_titles[[i]], title, sep="_")
    signature <- cell_signatures[[i]]
    x <- make_human_heatmap(exp_matrix, signature, title, average = FALSE)
    # dev.off()
  })
}
# Run function for each signature list (before, after ANOVA, and after ANOVA+Tukey)
make_cell_signature_heatmap(significant_hits, "significant_hits") #unclustered #cluster_cols = F






