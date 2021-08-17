rm(list=ls())

## 2. Ridgeplots to boxplots 
kazutils::load_packages(c("ggpubr", "ggplot2", "gtools", "plyr"))

# Gene sets 
sets <- c("HR", "NHEJ")

# # Meta data
ES_df <- readRDS("2a_enrichment_metadata.rds")
# unique(ES_df$phase) is c("S/G2/M", "M/G1/S")

# ANOVA functions
source("../functions_ANOVA_Tukey.R")

# Plot parameters
lvl.colors <-  c("#D55E00", "#0072B2")
font_size <- 15
line_size <- 1
pval.test <- "t.test"
pval.label <- "p.format"

# Plot for each set
for (set in sets){
  # Data frame
  df <- ES_df %>%
    mutate(set = ES_df[,set]) %>%
    .[, c("phase", "broadcelltype", "set")]
  df$broadcelltype <- plyr::mapvalues(df$broadcelltype, 
                                      from=c("Luminal progenitor", "Basal", "Mature luminal"), 
                                      to=c("LP","BC","ML"))
  # Column names
  colnames(df) <- c("variable", "level", "value")
  phases <- levels(df$variable)
  
  # ANOVA stats
  stats_anova <- lapply(phases, function(x){
    get_ANOVA_Tukey_text(df[df$variable ==x, 2:3], x)
  })
  
  # Make combinations of cell types for stats_compare_means
  # Make list of unique elements
  ele <- unique(as.character(df$level))
  # Make list of combinations (order doesn't matter) for p-values
  comb <- gtools::combinations(n = length(ele), r = 2, v = ele, repeats.allowed = F) %>% # gtools
    split(., seq(nrow(.)))
  # Star height relative to bars
  vjust <- ifelse(pval.label == "p.signif", 0.5, -0.1)
  
  # Plot
  p <- ggplot(df, aes(x = level, y = value, fill = level, color = level)) +
    # scale_x_reordered()+
    facet_wrap(~variable,  scales = "free_x")+
    labs(
      title = sprintf("%s enrichment across Revelio phases in MEC single cell", set),
      y = "escape::enrichIt score",
      x = "",
      subtitle = sprintf("stats: above brackets = %s, in caption = ANOVA/tukey; violin line = mean", pval.test),
      caption = paste(unlist(stats_anova), collapse = "\n")
    ) +
    theme_classic()+
    theme(
      #plot.subtitle = element_text(size=font_size),
      # Ticks
      axis.ticks = element_line(colour = "black", size = line_size), # increase the tick thickness)
      axis.ticks.length = unit(.25, "cm"),
      # Axes labels
      axis.text = element_text(colour = "black", size = font_size),
      axis.text.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0), hjust = 1, vjust = 1, angle = 45), # increase space between x axis title and labels
      axis.text.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
      # axes tick labels
      axis.title = element_text(colour = "black", size = font_size, face = "bold"), # axes title labels
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), # increase space between x axis title and labels
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
      legend.position = "none"
    ) +
    stat_compare_means(
      method = pval.test, comparisons = comb, na.rm = T, vjust = vjust, hide.ns = T,
      label = pval.label, size = font_size / 4, bracket.size = 0.1
    )
  p1 <- p + geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(alpha=0.5, size =0.5,width = 0.2)
  p2 <- p + 
    geom_violin(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(alpha=0.5, size =0.5,width = 0.2)+
    stat_summary(fun.y=mean, geom="crossbar", size=0.1, color="black")
  # Plot and save to file
  pdf(sprintf("%s/2c_boxplot_%s.pdf", out_dir, set), width = 7)
  print(p1)
  print(p2)
  dev.off()
}

