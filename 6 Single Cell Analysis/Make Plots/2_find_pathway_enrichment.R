#' Purpose:
#' To calculate enrichment of set of damage genes in each cell using escape package.

# Remove all objects in environment except b2 (b2 = Seurat object)
rm(list=ls()[!ls()%in%"b2"])

# Libraries
library(kazutils)
load_packages(c("Seurat", "dplyr", "ggplot2", "escape", "GSEABase", "dplyr", "cowplot"))
source("../plot_ridgeEnrichment.R")

# Out directory
out_dir <- kazutils::create_folder("Plots")

# # Data 
# b2 <- readRDS("../annotated_Seurat_b2.rds")

# Gene lists 
cell_colors <- readRDS("../colors_broadcelltype.rds")
Revelio_colors <- readRDS("../colors_Revelio_ccPhase.rds")
load("../genes_HK_NHEJ_HR.RData")
genes_Rad51 <- c("RAD51B", "RAD51C", "RAD51D") # from rownames(b2)[grep("RAD51", rownames(b2))]

# Create GeneSet object from the gene lists 
sets <- c("RAD51_set", "HR", "NHEJ")
GS_RAD51 <- GSEABase::GeneSet(genes_Rad51, setIdentifier = "RAD51_set", setName = "RAD51_set")
GS_HR <- GeneSet(genes_HR, setIdentifier = "HR", setName = "HR")
GS_NHEJ <- GeneSet(genes_NHEJ, setIdentifier = "NHEJ", setName = "NHEJ")

# Create a collection of GeneSet objects
GSC_damage <- GSEABase::GeneSetCollection(GS_RAD51, GS_HR, GS_NHEJ)

### Subset b2 to only S, G2.M, G2 (proliferating = p) -------------------------
b2$Revelio_ccPhase %>% unique
b2_p <- subset(b2, Revelio_ccPhase %in% c("S","G2", "G2.M"))
dim(b2_p) # 26027 genes 2893 cells

# Enrich for the signatures on the raw counts for each group of genesets
ES <- escape::enrichIt(obj = b2, gene.sets = GSC_damage, groups = 1000) 

# Add enrichment scores to seurat object metadata 
b2 <- Seurat::AddMetaData(b2, ES)

# Save Data 
b2 <- saveRDS("../annotated_Seurat_b2.rds")