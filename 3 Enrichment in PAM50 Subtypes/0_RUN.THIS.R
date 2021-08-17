# Nov 12, 2019
# For HK, look for enrichment of total proteome signatures (FC > 3) in METABRIC cohort
FC <- 5
# Just run these three lines to make the 3 box plots and heat map
source("1_compute.ssGSEA.R") # correlate human signatures to METABRIC Breast Cancer subtypes
source("2_make.violin.and.box.plots.R") # show enrichment of each signature (Bc, LM, LP) in different BC subtypes 
source("3_make.heatmap.R") # visualize all ssGSEA enrichment scores across the BC METABRIC cohort

# Credit: original code was from Luis/Miquel, just took out survival analysis