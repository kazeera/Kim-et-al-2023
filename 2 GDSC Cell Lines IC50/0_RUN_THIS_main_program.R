#' Purpose: 
#' To perform interrogation of GDSC data. Look at the enrichment of our signatures in 
#' ..breast cancer cell lines and correlate to IC50 values.
#' 
#' Signature files were created in folder 1.

## Load libraries
library(openxlsx)
library(GSVA) # to compute ssGSEA
library(scales) # for rescaling scores
library(pheatmap) # to make heat map
library(RColorBrewer) # to pick heat map colors/pallettes
library(openxlsx) # to read annotation files from Excel
library(foreach) # to do drug correlations
library(reshape2)
library(scales)

## User defined analysis parameters:
# For human, change also in script 4a
tissue_type <- "breast" # can be "ovary", etc.
label <- "human_signatures"
signature_filename <- "human.signatures_FC5.pval0.05_premeno.only.RData"

# For mouse
# tissue_type <- "breast"
# label <- "mouse_signatures"
# signature_filename <- "hm.homolog.EP.mm.mamm.sig_FC3_pval0.05.RData"

TOP_WHAT <- 10 # this parameter defines how many cell lines that were enriched for BC AND LP 
               # ex. if TOP_WHAT = 5, we only look at the top 5 cell lines from basal and top 5 from LP 
               # define as 0 if you want all cell lines to be included

# Required files (A,B,C):
# A. IC 50 of drugs
load( "IC50_breast.cell.lines_HK.drugs.RData") # made in part 3
# B. ssGSEA of signature enrichment in cell lines
load(sprintf("ssGSEA_values_%s_cell_line_ssGSEA_%s_unscaled.RData", tissue_type, label)) # made in part 2a
# C. signatures is a list of 3 elements- BC, LM, LP - each element contains a vector with gene names
load(signature_filename)

## Pipeline:
# Create output folder
output_folder <- sprintf("Human Breast Cell Lines IC50 (FC%s, p%s, %s) Top%s", FC, pval_thres, label, TOP_WHAT)
dir.create(output_folder)

# Note: scripts that are commented out have data files that are not included/needed in this folder or not used for main analysis.
# source("1a.prepareExpressionMatrix.R")
load("Cell_line_RMA_proc_basalExp.RData") # made in part 1.a.
# source("1b.prepareSignatures.R") # only if files are in xlsx format   
source("2a.computeGSVA.R")     
source("2b.make_ssGSEA_cell.line_sig_heatmap.R") # Main Analysis for Hyeyeon - produces the BC LP heat map figure and excel file
# source("3.prepareIC50_mamm.R")                                                 
source("4a.computePCC.R")  
# source("4b.do.wilcox.test.R") - not used 

# Clear R environment
rm(list=ls()) 
