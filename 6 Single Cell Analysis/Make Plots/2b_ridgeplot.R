rm(list=ls())

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

# Meta data
ES_df <- readRDS("2a_enrichment_metadata.rds")
# unique(ES_df$phase) is c("S/G2/M", "M/G1/S")
# https://cran.r-project.org/web/packages/ggridges/vignettes/gallery.html

sets <- c("RAD51_set", "HR", "NHEJ")
plots <- lapply(sets, function(set){
  # Data frame
  df <- ES_df %>%
    mutate(set = ES_df[,set]) %>%
    .[, c("set", "broadcelltype", "phase")]
  
  # Plot
  ggplot(df, aes(x = set, y = broadcelltype, fill = phase, color = phase)) +  
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
  ggsave(sprintf("%s/2b_ccPhase_ridgeplot_%s.pdf", out_dir, set), height = 3, width = 6)#,height = 10)
  
  # # Return
  return(p)
})

# Data frame
df <- ES_df %>%
  mutate(set = ES_df[,'NHEJ']) %>%
  .[, c("set", "broadcelltype", "phase")]

source("../functions_ANOVA_Tukey.R")


df <- ES_df %>%
  mutate(set = ES_df[,'HR']) %>%
  .[, c("set", "broadcelltype", "phase")]
paste(get_ANOVA_Tukey_text(df[df$broadcelltype == "LP",c("phase", "set")]))
t.test()