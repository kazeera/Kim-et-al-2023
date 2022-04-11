#' Purpose: 
#' Create cell-specific proteomic signatures by looking at protein expressions and 
#' ..subset based on significance and fold change thresholds defined in script 0.RUN_THIS
#' 
#' Note: July 5, 2019 ---> revised Nov 15, 2019 for HK's work 
#' This script was made for Mathepan's metabolic proteome analysis but we can use it to make total proteome signatures
#' 
#' Revisions on Nov 15, 2019:
#' Added lines so we only look at premeno samples
rm(list=ls())

# Import required library
library(openxlsx)

# Load functions and data
load("h.pint.pheno_combat.proteome_premeno.only.RData")
source("makeSignatures_functions.R")

# Since old code was tailored to all samples, reassign objects without post.meno samples to old variable names
h.pint.combat <- h.pint.combat_matrix_no.post.meno
pint.pheno <- pint.pheno_no.meno

# Analysis ----
# Prepare expression matrices
exp_matrix <- h.pint.combat
# average
h.pint.combat.avg <- average_hmpint_columns(h.pint.combat) #make matrix with average
exp_matrix_avg <- h.pint.combat.avg

# Run ANOVA and Tukey's test on all proteins
# Make df with gene name
significance_df <- compute_ANOVA_tukey_pvals(exp_matrix)
rownames(significance_df) <- rownames(h.pint.combat)

# 2. Calculate log 2 FC of each cell type vs the other
log2FC_df <- compute_logFC_3groups(exp_matrix_avg)
rownames(log2FC_df)<- rownames(h.pint.combat)

# Get data frames for seperate analyses
subset_dfs <- lapply(c("BC", "LM", "LP"), function(ID){
  subset_df(significance_df, log2FC_df, ID)
})

# Save workspace as RData
save.image(file = "human_BC.LM.LP_combat.pint_significance.log2FC.tables_premeno.only.RData")
