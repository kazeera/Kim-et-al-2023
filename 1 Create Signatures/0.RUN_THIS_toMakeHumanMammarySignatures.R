rm(list=ls())
# Load data 
load("human_BC.LM.LP_combat.pint_significance.log2FC.tables_premeno.only.RData")

# 1. Compute significance and fold-change values between BC, ML, LP
source("1_makeSignificance_and_log2FoldChange_tables.R")

# 2a. Define thresholds
pval_thres <- 0.05
FC <- 3
log2FC_thres_up <- log-2(FC)

# 2b. Get proteins that meet thresholds and create signature lists 
source("2_defineThresholds_and_makeSignatures.R")
