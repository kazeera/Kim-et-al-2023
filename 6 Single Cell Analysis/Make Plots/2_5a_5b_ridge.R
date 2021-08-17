# This script subsets 
# - Revelio cell cycle phase to proliferating (p) and non proliferating (np)
# - finds enrichment of damage genes
# - then compares the two in scripts labelled 5b, 5c

# Libraries
library(kazutils)
load_packages(c("Seurat", "dplyr", "ggridges", "ggplot2", "escape", "GSEABase", "dplyr", "cowplot"))
source("../plot_ridgeEnrichment.R")

# Out directory
out_dir <- kazutils::create_folder("Final Plots")

# Data 
b2 <- readRDS("../annotated_Seurat_b2.rds")

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
ES_p <- escape::enrichIt(obj = b2_p, gene.sets = GSC_damage, groups = 1000) 

# Add enrichment scores to seurat object metadata 
b2_p <- Seurat::AddMetaData(b2_p, ES_p)

# Subset Seurat metadata to just plot certain clusters on ridgeplot below
ES_p <- data.frame(b2_p[[]], Idents(b2_p))

### Subset b2 to only M.G1, G1.S (non-proliferating = np) -------------------------
b2$Revelio_ccPhase %>% unique
b2_np <- subset(b2, Revelio_ccPhase %in% c("M.G1", "G1.S"))
dim(b2_np) # 26027 genes 2893 cells

# Enrich for the signatures on the raw counts for each group of genesets
ES_np <- escape::enrichIt(obj = b2_np, gene.sets = GSC_damage, groups = 1000) 

# Add enrichment scores to seurat object metadata 
b2_np <- Seurat::AddMetaData(b2_np, ES_np)

# Subset Seurat metadata to just plot certain clusters on ridgeplot below
ES_np <- data.frame(b2_np[[]], Idents(b2_np))

### Combine n and np
ES_df <- c("np", "p") %>%
  lapply(function(x){
    # Use proper ES
    if(x=="np")
      ES <- ES_np
    else
      ES <- ES_p
    # Rename columns
    cbind(ES[,c("broadcelltype", "Revelio_ccPhase", sets)], phase = x)
  })
ES_df <- do.call(rbind.data.frame, ES_df)
ES_df$broadcelltype <- as.character(ES_df$broadcelltype)
ES_df$phase <- ES_df$phase %>%
  as.character %>% 
  plyr::mapvalues(from=c("p","np"), to=c("S/G2/M", "M/G1/S")) %>%
  factor(levels = c("S/G2/M", "M/G1/S"))
ES_df$Revelio_ccPhase <- as.character(ES_df$Revelio_ccPhase)

saveRDS(ES_df, "5a_enrichment_metadata.rds")

# # Meta data
ES_df <- readRDS("5a_enrichment_metadata.rds")
# unique(ES_df$phase) is c("S/G2/M", "M/G1/S")
# https://cran.r-project.org/web/packages/ggridges/vignettes/gallery.html

plots <- lapply(sets, function(set){
  # Data frame
  df <- ES_df %>%
    mutate(set = ES_df[,set]) %>%
    .[, c("set", "broadcelltype", "phase")]
  
  # Plot
  p <- ggplot(df, aes(x = set, y = broadcelltype, fill = phase, color = phase)) +  
    geom_density_ridges(
      jittered_points = TRUE, scale = .95, rel_min_height = .01,
      point_shape = "|", point_size = 3, 
      position = position_points_jitter(height = 0),
      quantile_lines = T, quantile_fun = function(x,...)mean(x), size = 1
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("#D55E0050", "#0072B250")) +
    scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
    scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(
      override.aes = list(
        fill = c("#D55E00A0", "#0072B2A0"),
        color = NA, point_color = NA)
    )) +
    labs(
      title = set,
      x = "score") +
    theme_ridges(center = TRUE) +
    theme(legend.position = "bottom")
  # Save plot
  # ggsave(sprintf("%s/5b_MEC_ridgeplot_2sets_%s.jpeg", out_dir, set), height = 3, width = 6)#,height = 10)
  ggsave(sprintf("%s/5b_MEC_ridgeplot_2sets_%s.pdf", out_dir, set), height = 3, width = 6)#,height = 10)
  
  # # Return
  return(p)
})