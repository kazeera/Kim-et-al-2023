#' Purpose:
#' To calculate ssGSEA scores to see enrichment of gene lists (ie signatures) in GDSC RNA expression data
#' 
#' Credit: Pujana lab for conception, Luis for original code

# Load R libraries
library(GSVA) # to compute ssGSEA
# library(scales) # for rescaling scores
library(pheatmap) # to make heat map
library(RColorBrewer) # to pick heat map colors/palettes
library(openxlsx) # to read annotation files from Excel

# Load GDSC data  
load("Cell_line_RMA_proc_basalExp_clean.RData")

# Load signature data
# # load(list.files(pattern="signatures|RData"))
# load(signature_filename)

# Get cell line names from annotations table
cell_line_ann <- read.xlsx("TableS1E Cell Lines clean.xlsx")

# Subset to only cell lines of interest (specific to tissue type)
cl_breast_indices <- which(cell_line_ann$GDSC.Tissue.descriptor.2 == tissue_type) #tissue_type defined in script 0
cell_lines_of_interest <- cell_line_ann$COSMIC.identifier[cl_breast_indices]

# number of cell lines that are in cell lines of interest
cols_to_keep <- which(colnames(exp_matrix) %in% cell_lines_of_interest)
length(cols_to_keep)

# Subset exp_matrix to these columns only
exp_matrix_br <- exp_matrix[,cols_to_keep]
# save(exp_matrix_br, file=sprintf("Cell_line_RMA_proc_basalExp_clean_%s.RData", tissue_type)) 


# Compute ssGSEA
ssgseaValues = gsva(as.matrix(exp_matrix_br), 
                    signatures, 
                    method='ssgsea', 
                    kcdf = 'Gaussian', 
                    verbose = FALSE)

# library(scales)
# plot_matrix <- apply(t(ssgseaValues_pdx_raw), 2, function(x) rescale(x, c(-3,3)))

# Make heat map of enrichment scores
pheatmap(t(ssgseaValues),
         scale = "none",
         cellwidth = 50,
         cellheight = 10,
         # color = brewer.pal(9, "RdYlGn"),
         filename = sprintf("%s/heatmap_%s_cell_line_ssGSEA_%s%s_unscaled_premeno.png", output_folder, tissue_type, label, FC))
save(ssgseaValues, file = sprintf("%s/ssGSEA_values_%s_cell_line_ssGSEA_%s_unscaled.RData", output_folder, tissue_type, label))

save.image(sprintf("%s/2.ssGSEA.RData", output_folder))

