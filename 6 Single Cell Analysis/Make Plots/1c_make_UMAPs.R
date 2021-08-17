# Cell count
# table(b2$broadcelltype)
# Mature luminal Luminal progenitor              Basal 
# 2172               3364                         1272 

rm(list=ls()[!ls()%in%"b2"])
# Final plots
library(kazutils)
load_packages(c("Seurat"))

# # Data 
b2 <- readRDS("../annotated_Seurat_b2.rds")

# Out directory
out_dir <- kazutils::create_folder("Final Plots")

# Colors
cell_colors <- readRDS("../colors_broadcelltype.rds")
names(cell_colors) <- c("BC", "LP", "ML")
Revelio_colors <- readRDS("../colors_Revelio_ccPhase.rds")
Phase_colors <- readRDS("../colors_Revelio_ccPhase.rds")

# New identities
Idents(b2) <- b2$broadcelltype2

## 1. UMAPs
p1 <- DimPlot(b2, group.by = "Revelio_ccPhase", cols = Phase_colors, reduction = "umapharmony", pt.size = 0.1)+
  theme(legend.position = "bottom")
p2 <- DimPlot(b2, group.by = "broadcelltype2", cols = cell_colors, reduction = "umapharmony", pt.size = 0.1)+
  theme(legend.position = "bottom")

# Plot UMAP of all phases
pdf(sprintf("%s/1c_umap_total_ccPhase.pdf", out_dir))
p1
p2
dev.off()

# Plot UMAP of all phases
pdf(sprintf("%s/1c_umap_total_ccPhase_noLegend.pdf", out_dir))
p1 + NoLegend()
p2 + NoLegend()
dev.off()


# Plot UMAPs
out_dir <- create_folder("Final Plots")
pdf(sprintf("%s/1c_umap_total_GMNN.pdf", out_dir))
FeaturePlot(b2, "GMNN", reduction = "umapharmony", col = c("gray18", "yellow"), pt.size = 0.5) +DarkTheme()
FeaturePlot(b2, "GMNN", reduction = "umapharmony", pt.size = 0.5)
dev.off()
