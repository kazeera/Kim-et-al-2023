rm(list=ls())
load("human_BC.LM.LP_combat.pint_significance.log2FC.tables_premeno.only.RData")

# source("1_makeSignificance_and_log2FoldChange_tables.R")
# Assign cell type to row "significant hits"
# define thresholds
pval_thres <- 0.05
FC <- 3
log2FC_thres_up <- log2(FC)
source("2_defineThresholds_and_makeSignatures.R")
