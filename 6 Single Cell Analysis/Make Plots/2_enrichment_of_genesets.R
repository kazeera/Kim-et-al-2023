# This script subsets 
# - Revelio cell cycle phase to proliferating (p) and non proliferating (np)
# - finds enrichment of damage genes
# - then compares the two in scripts labelled 5b, 5c
rm(list=ls()[!ls()%in%"b2"])

# Libraries
library(kazutils)
load_packages(c("Seurat", "dplyr", "ggridges", "ggplot2", "escape", "GSEABase", "dplyr", "cowplot"))
source("../plot_ridgeEnrichment.R")

# Out directory
out_dir <- kazutils::create_folder("Final Plots")

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

# Enrich for the signatures on the raw counts for each group of genesets
ES <- escape::enrichIt(obj = b2, gene.sets = GSC_damage, groups = 1000) 

# Add enrichment scores to seurat object metadata 
b2 <- Seurat::AddMetaData(b2, ES)
all(ES$HR == b2$HR) # TRUE

# Data 
saveRDS(b2, "../annotated_Seurat_b2.rds")

# Subset Seurat metadata to just plot certain clusters on ridgeplot below
ES_df <- data.frame(b2[[]], Idents(b2))

ES_df$broadcelltype <- ES_df$broadcelltype %>%
  as.character %>% 
  plyr::mapvalues(from=c("Mature luminal","Luminal progenitor","Basal"), to=c("ML", "LP", "BC")) %>%
  factor(levels = c("LP", "ML", "BC"))

ES_df$phase <- ES_df$Revelio_ccPhase %>%
  as.character %>% 
  # M.G1, G1.S (non-proliferating) = "M/G1/S"
  # S, G2.M, G2 (proliferating) = "S/G2/M"
  plyr::mapvalues(from=c("M.G1", "G1.S","S", "G2", "G2.M"), to=c("M/G1/S","M/G1/S","S/G2/M","S/G2/M","S/G2/M")) %>%
  factor(levels = c("S/G2/M", "M/G1/S"))

# > table(ES_df$phase)
# S/G2/M M/G1/S 
# 2893  + 1914 = 4807
# NA 2001

ES_df <- ES_df[!is.na(ES_df$phase), c("broadcelltype", "phase", "RAD51_set", "HR", "NHEJ")]

saveRDS(ES_df, "2a_enrichment_metadata.rds")
