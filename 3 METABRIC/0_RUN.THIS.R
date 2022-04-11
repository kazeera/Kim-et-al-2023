#' Purpose:
#' To interrogate METABRIC Breast Cancer (BC) data by running ssGSEA that finds enrichment of signatures (gene sets)
#' ..in BC subtypes (Pam-50 + Claudin Low). RNA expression data is from the METABRIC study from RData file in the folder. 

# For HK, look for enrichment of total proteome signatures in METABRIC cohort
FC <- 5

# Just run these three lines to make the 3 box plots and heat map
source("1_compute.ssGSEA.R") # correlate human signatures to METABRIC Breast Cancer subtypes
source("2_make.violin.and.box.plots.R") # show enrichment of each signature (Bc, LM, LP) in different BC subtypes 
source("3_make.heatmap.R") # visualize all ssGSEA enrichment scores across the BC METABRIC cohort

#' Credit: Pujana lab for conception, adopted from Luis' scripts - modified so that relevant code in 3 scripts compiled into 1 
# Nov 12, 2019