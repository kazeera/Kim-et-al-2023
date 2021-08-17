plot_ridgeEnrichment <- function(ES, gene_set_id, group_by, colors, plot_title=""){
  library(ggridges)
  library(escape)
  p <- escape::ridgeEnrichment(ES, gene.set = gene_set_id, group = group_by, add.rug = TRUE,
                  colors = colors) +
    theme_bw()+
    # theme(axis.text.x= element_blank(),
    #       plot.title = element_text(lineheight=.8, face="bold"),
    #       axis.text = element_text(size = 35, color = "black", family = "sans"),
    #       axis.title.x=element_blank(),
    #       axis.line = element_line(colour = 'black'),
    #       axis.ticks = element_line(colour = 'black'),
    #       axis.title = element_text(size = 15, face = 'bold'),
    #       axis.title.y = element_blank(),
    #       legend.title=element_text(size=10,face = 'bold'),
    #       legend.text=element_text(size=10),
    #       panel.background = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.grid.major = element_line(colour = 'white'),
    #       plot.background = element_blank(),
    #       legend.key = element_rect(color = "transparent", fill = "transparent")) +
    ggtitle(plot_title)+
    geom_density_ridges(quantile_lines = T, quantile_fun = function(x,...)mean(x), size = 1.25) #ggridges
  print(p)
}
