# Nov 8, 2019
rm(list=ls())
# This script subsets the mouse proteome based on a genes of interest list

# Input: HK's gene list and mouse proteome data
#' **** HK wanted to only pull out DDR pathway genes in the proteome using R, 
#'      so she made expression for the other genes "NaN"
#'      I use this to make a list of genes to use in the prepared expression matrix
#'      


## Call libraries
library(cluster) #for diana clustering
library(scales)
library(pheatmap)
library(RColorBrewer)

source("plot_m_heatmap_functions.R")

# Load mouse proteome data from step one. We only need z.scores
load("mouse_proteome.RData")
# Transpose matrix
m.z.scores <- t(m.z.scores)

# # Reassign row names as 6-character PROBE ID
# rownames(m.z.scores) <- substring(as.character(m.ibaq.original$X),1,6)

# Make new column with 6-character PROBE ID
m.z.scores <- as.data.frame(m.z.scores)
m.z.scores$PROBE.ID <- substring(as.character(m.ibaq.original$X),1,6)

# Remove all envt variables except mouse z scores
rm(list=ls()[-which(ls()=="m.z.scores")])
head(m.z.scores)

## GET GENES OF INTEREST ---------
# Input: LFQ-adjusted IBAQ values (total proteome)## Read in proteome data
cytoscape_original <- read.delim(file = "dna_repair_genes_cytoscape.txt", header = T)
# > head(ddr.cytoscape_original[,c(19,28:30)])
# EM17_GS_DESCR EM17_pvalue..nc.                                                                   name selected
# 1        Dual incision in TC-NER               NA        DUAL INCISION IN TC-NER%REACTOME DATABASE ID RELEASE 63%6782135     TRUE
# 2                DNA replication               NA                                        DNA REPLICATION%GOBP%GO:0006260     TRUE

cytoscape_ddr <- cytoscape_original[which(cytoscape_original$selected),]

# Extract the gene names 
# First split the probe ids with "|"
ddr_genes <- strsplit(as.character(cytoscape_ddr$EM17_Genes), split = "\\|")
names(ddr_genes) <- gsub("%.*", "", cytoscape_ddr$name)


# Then make them into a vector and remove duplicates
genes_of_interest <- ddr_genes$`DNA DOUBLE-STRAND BREAK REPAIR`
genes_of_interest <- unlist(ddr_genes, use.names = F)
genes_of_interest <- unique(genes_of_interest)


DSB_df <- plot_genes_of_interest(ddr_genes$`DNA DOUBLE-STRAND BREAK REPAIR`, "DSB")
DDR.cell.resp <- plot_genes_of_interest(ddr_genes$`CELLULAR RESPONSE TO DNA DAMAGE STIMULUS`, "cell.response.to.DDR")



library(openxlsx)

wb <- createWorkbook()
write_df_to_Excel(current_sheet = "DSB exp", my_df = DSB_df, wb, incl_rowNames = TRUE)
write_df_to_Excel(current_sheet = "cell response to DDR exp", my_df = DDR.cell.resp, wb, incl_rowNames = TRUE)
write_df_to_Excel(current_sheet = "Total mouse proteome", my_df = m.z.scores, wb, incl_rowNames = TRUE)
write_df_to_Excel(current_sheet = "Cytoscape DNA Damage Cluster", my_df = cytoscape_ddr, wb, incl_rowNames =  FALSE)
write_df_to_Excel(current_sheet = "Cytoscape Original", my_df = cytoscape_original, wb, incl_rowNames = FALSE)
excel_filename <- sprintf("%s_HK_DDR.cytoscape.genes_mouse.mamm.proteome.z.scores.xlsx", format(Sys.Date(), "%y%m%d"))
saveWorkbook(wb, file = excel_filename,overwrite = T)

save.image("workspace.RData")
q()
