#' Purpose: 
#' To read in ibaq-adj lfq proteome expression data, and make the following modifications:
#' a) make all gene unique in the matrix, 
#' b) impute zeros with a random value between 1 and 1.5 (for later stats purposes)
#' c) log2 the data - transform data (log2) so it will be normally distributed.
#' Save backups as files and RData
#' 
#' Credits: modified Alison Casey's code. 

# Clear R environment
rm(list=ls())

# Import R libraries
library(cluster) #for diana clustering
library(scales)
library(pheatmap)
library(RColorBrewer)

# Import functions
source("plot_m_heatmap_functions.R") # plot mouse (m) heatmap

# Input: LFQ-adjusted IBAQ values (total proteome)## Read in proteome data
m.pint <- read.csv(file = "casey_IBAQadjustedLFQ_mouse_proteome_sort.csv")

# Check for the following.. empty gene names, duplicate gene names, expression value of 0 (non-imputed)
any(m.pint$gene == "") #TRUE
any(duplicated(m.pint$gene) & m.pint$gene != "") #TRUE
any(m.pint == 1) #TRUE - raw data, alison imputed 0s with 1s

# Make gene names unique, e.g. 3 instances of PKM will become PKM, PKM.1 and PKM.2
m.pint$gene <- make.unique(as.character(m.pint$gene), sep=".")
any(duplicated(m.pint$gene)) #FALSE

# Sort matrix based on gene symbol (alphabetically)
m.pint <- m.pint[order(m.pint$gene),] 
m.pint -> m.ibaq.original

# Define helper functions to create Z score matrix
# This function replaces all 1s in a vector x with a random value between start_index and last_index
replace_ones_with_double <- function(x, start_index, last_index){
  x[x == 1]<- runif(sum(x == 1), start_index, last_index)
  return(x)
}

# This function computes and returns matrix with z scores for heatmap plotting
compute_z_scores <- function(matrix){
  # Formula for z score: 
  #         [element x - mean of row x is in]/ standard deviation of row x is in
  return (apply(matrix, 1, function(x) (x - mean(x)) / sd(x)))
}

# Create new data frame with relevant sample information
m_exp_matrix <- data.frame(m.pint$e.b.1, m.pint$e.b.2, #B E
                           rowMeans(m.pint[c('e.p.b.1', 'e.p.b.1.rep')], na.rm=TRUE), rowMeans(m.pint[c('e.p.b.2', 'e.p.b.2.rep')], na.rm=TRUE), #B EP #take average of replicates
                           m.pint$e.ap.1, m.pint$e.ap.2, #LP E
                           m.pint$e.p.ap.1, m.pint$e.p.ap.2, #LP EP
                           m.pint$e.pr.l.1, m.pint$e.pr.l.2, #LM E
                           m.pint$e.p.pr.l.1, m.pint$e.p.pr.l.2) #LM EP
# Rename rows of matrix as genes
rownames(m_exp_matrix) <- m.pint$gene


# Define columns annotations (sample descriptions)
cell_type_samples <- rep(c("B", "LP", "LM"), each = 4)
hormone_samples <- rep(c("E","EP"), each = 2, times = 3)
replicate <- rep(1:2, times = 6)

# Rename columns with sample information (celltype and hormone treatment)
colnames(m_exp_matrix) <- paste(cell_type_samples, hormone_samples, replicate, sep="_")

m_exp_matrix -> m.ibaq

#impute 1s with random number between 1 and 1.5 
m_exp_matrix <- replace_ones_with_double(m_exp_matrix, 1, 1.5)
m_exp_matrix -> m.imputed
# log2 all values to make normal distribution
m_exp_matrix <- log2(m_exp_matrix)
m_exp_matrix -> m.imputed.log2
#convert values into z scores
m_exp_matrix_z <- compute_z_scores(m_exp_matrix)
m_exp_matrix_z -> m.z.scores

# Save R data file
save(m.ibaq.original, m.ibaq, m.imputed, m.imputed.log2, m.z.scores, file = "mouse_proteome.RData")

# Make dendograms -by Alison Casey
## a) For samples (columns)
dist_samples <- as.dist(1-cor(m_exp_matrix,
                              #use = "pairwise.complete.obs",
                              method = "pearson"))
hclust_samples <- hclust(dist_samples)
# Plot
png("sample_dendogram.png")
plot(as.dendrogram(hclust_samples)) 
dev.off()

## b) For proteins (rows)
dist_proteins <- as.dist(1-cor(t(m_exp_matrix),
                               #use = "pairwise.complete.obs",
                               method = "pearson"))
# be prepared - this next line of code takes like 5 minutes
hclust_proteins <- hclust(dist_proteins) #DIvisive ANAlysis Clustering

# Save RData
save(hclust_samples, dist_samples, hclust_proteins, dist_proteins, file = "mouse_dendograms.RData")

# Plot heatmap
plot_mouse_heatmap(m_exp_matrix_z, "mouse_heatmap.png", NA, cell_type_samples, hormone_samples,hclust_proteins, hclust_samples)

# Save final matrix with only EP 
save(m.imputed.log2, file = "mouse_imputed.log2.RData")
