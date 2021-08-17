# Helper functions to sort facet by median
# https://stackoverflow.com/questions/47090344/how-to-properly-sort-facet-boxplots-by-median/47091845
reorder_within <- function(x, by, within, fun = median, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

# df2 in this order: 1) box 2) facet by 3) numeric value
# cols = named vector of colors for boxes
facet_boxplot <- function(df2, cols = NA, plot_title = "",add_mean = T, sub_title = "RM", x_lab = "Phase", y_lab = "Expression"){
  colnames(df2) <- c("var", "group", "value")
  p <- ggplot(df2, aes(x = reorder_within(var, value, group, median), y = value)) + 
    geom_boxplot(aes(fill = var), outlier.shape = NA, width = 5) + 
    geom_jitter(size = 3, width = 0.1, alpha = 0.7, pch = 21) +
    scale_x_reordered()+
    facet_wrap(~group,  scales = "free_x")+
    labs(title = plot_title,
         subtitle = sub_title, 
         x = x_lab,
         y = y_lab)+
    theme(panel.background = element_blank(),#remove background color and lines
          axis.line = element_line(colour = 'black', size = line_size), # increase the axis-line thickness and change the color to blac
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), #increase space between x axis title and labels
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          legend.title = element_blank())
  if(!is.na(cols)){
    p <- p +
      scale_color_manual(values = Seurat_phase_colors) +
      scale_fill_manual(values = Seurat_phase_colors)
  }
  if(isTRUE(add_mean)){
    p <- p+ 
      stat_summary(# Add a line for mean
        fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
        position = position_dodge(width = 0.3), width = 0.25, color = "black", size = 2
      )
  }
  print(p)
}

boxplot1 <- function(df, add_box = F, cols = NA, line_size = 1, plot_title = "",add_mean = T, sub_title = "RM", x_lab = "Phase", y_lab = "Scaled Expression"){
  colnames(df) <- c("var", "value")
  
  p <- ggplot(df, aes(x = var,  y = value, color=var)) 
  
  if(add_box){
    p <- p +
      geom_boxplot(size=1, na.rm = T, outlier.colour = NA) +
      geom_jitter(size = 2, width = 0.2, color = "black", fill = "black",  alpha = 0.7, pch = 21) 
  } else {
    p <- p +
      geom_jitter(aes(fill=var), size = 3, width = 0.2, alpha = 0.7, pch = 21) 
  }
  
  if(isFALSE(is.na(cols))){
    p <- p +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols)
  }
  
  p <- p+
    theme_light() +
    labs(    
      title = plot_title,
      subtitle = sub_title, 
      x = x_lab,
      y = y_lab)+ 
    stat_summary(# Add a line for mean
      fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
      position = position_dodge(width = 0.3), width = 0.25, color = "black", size = 2
    )
  
  print(p)
}
