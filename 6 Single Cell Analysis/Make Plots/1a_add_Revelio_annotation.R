#' Purpose:
#' To look for expression of genes in different cell cycle (cc) phases using Revelio package.
#' This script adds Revelio (cc package) ccPhase, ccAngle, etc annotations to Seurat cells meta data
#' 
#' Credit: Curtis McCloskey created, normalized and interrogated reduction mammoplasty single cell data in Mathepan's paper (https://github.com/mcclo/Mahendralingam-et-al.-Nat-Metab)
#' modified Mar 2, 2021

# Import required R packages
library(kazutils)
load_packages(c("Seurat", "plyr", "Revelio", "viridis", "HGNChelper"))
# library(scClustViz) # BiocManager::install("scClustViz")

# Load Seurat object # this data file was too large to upload
load("../Data/b2_obj_only.RData")

# Make a column with shortened cell type names
b2$broadcelltype2 <- plyr::mapvalues(b2$broadcelltype, from=c("Mature luminal", "Luminal progenitor", "Basal"), to=c("ML","LP","BC"))

# Cell cycle phases
phases <- c("G1.S", "S", "G2.M", "M.G1", "G2")

# Load revelio genes (need to update)
data("revelioTestData_cyclicGenes")

# Make list of updated phase gene markers 
phase_genes <- lapply(phases, function(phase){
  g1s <- revelioTestData_cyclicGenes[,phase]
  
  # use HGNhelper to identify incorrect symbols
  testg1s <- HGNChelper::checkGeneSymbols(g1s)
  
  # extract corrected info and make new matrix
  corr <- testg1s$Suggested.Symbol
  marker <- b2[intersect(corr, rownames(b2)),]
  diff <- setdiff(corr, rownames(marker))
  remove <- as.character(diff)
  # return corrected info
  corr[! corr %in% remove]
})
# Rename list 
names(phase_genes) <- phases
# lapply(phase_genes, length) # varying lengths

# Convert a list with various length vectors to a data.frame - placeholders become NA
reveliocorr <- data.frame(lapply(phase_genes, "length<-", max(lengths(phase_genes))))
# all(reveliocorr2[,2] == reveliocorr[,2], na.rm = T) #TRUE

# Load matrix from seurat object
mat <- as.matrix(GetAssayData(b2, slot = "counts"))

# Create Revelio object
myData <- createRevelioObject(rawData = mat, cyclicGenes = reveliocorr)

# Rename batchIDs to match expected input (ie. batch ID (all 1 batch) instead of cell ID)
myData@cellInfo[,'batchID'] <- factor(rep("group1", length(myData@cellInfo[,'batchID'])))

# Perform Revelio workflow steps
#might have to look at ?getCellCyclePhaseAssignInformation for cutoffs since it doesnt assign to every cell
myData <- getCellCyclePhaseAssignInformation(dataList = myData)
myData <- getPCAData(dataList = myData, boolPlotResults = TRUE)
myData <- getOptimalRotation(dataList = myData, boolPlotResults = TRUE) # takes a min or 2

# Add to Seurat object metadata 
# Revelio columns/parameters to add: 
rev_cols <- c("ccPhase", "ccAngle", "ccRadius", "ccTime", "ccPercentage", "ccPercentageUniformlySpaced", "ccPositionIndex")
rev_df <- myData@cellInfo[, rev_cols]
# Add in a loop
for(cc in rev_cols){
  b2 <- AddMetaData(b2, rev_df[,cc,drop=F], col.name = paste("Revelio", cc, sep="_"))
}

# Print how many cells have labelled cell cycle phases
sprintf("%s cells have phase info out of %s total cells", sum(rownames(rev_df) %in% colnames(b2)), ncol(b2))
# "4807 cells have phase info out of 6808 total cells"

win.graph()
DimPlot(b2, group.by = "Revelio_ccPhase", reduction = "umapharmony", pt.size = 1)


# Reassign identities
Idents(b2) <- "broadcelltype"

b2$Revelio_ccPhase <- factor(b2$Revelio_ccPhase, levels = c("M.G1", "G1.S", "S", "G2", "G2.M"))

# Save RData 
revelio_ob <- myData
save(revelio_ob, file = sprintf("../1_RM_Revelio_ob.RData"))
saveRDS(b2, "../annotated_Seurat_b2.rds")

# Colors for points = cell cycle phase
# Revelio_colors <- RColorBrewer::brewer.pal(6, "Spectral")
# Revelio_colors <- c(M.G1="#8F1713", G1.S="#BA5511", S="#C9C126", G2="#3F782C", G2.M="#4871CF")
Revelio_colors <- viridis(length(levels(b2$Revelio_ccPhase)))
names(Revelio_colors) <- levels(b2$Revelio_ccPhase)
# scales::show_col(Revelio_colors)
saveRDS(Revelio_colors, "../colors_Revelio_ccPhase.rds")


