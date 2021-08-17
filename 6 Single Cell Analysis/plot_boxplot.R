#' Creates single box plot.
#'
#' @family plotting
#' @param df A data frame in long format with 2-3 columns: 1) Box or level = factor/character, 2) Value = numeric, 3) Dots/points (color, optional) = factor/character.
#' @param mean_line Should we draw a line on the boxes to show mean?
#' @param colors  A vector of colors to manually specify box colors.
#' @param xlab X axis label.
#' @param ylab Y axis label.
#' @param pval.test String corresponding to method parameter in \code{\link[ggpubr]{stat_compare_means}}. Allowed values are "t.test" and "wilcox.test".
#' @param pval.label String corresponding to label parameter in \code{\link[ggpubr]{stat_compare_means}}. Allowed values are "p.signif" (stars) and "p.format" (number).
plot_boxplot <- function(df,colors,
                         title = "", subtitle = "", caption ="", add_Tukey = T, 
                         pval.test = NULL, pval.label = "p.format",font_size = 15,
                         xlab = NULL, ylab = "Scaled Expression",
                         box = FALSE, mean_line = FALSE){
  library(dplyr)
  library(gtools)
  library(ggplot2)
  library(ggpubr)
  
  # x axis title
  xlab <- ifelse(is.null(xlab), colnames(df)[1], xlab)
  # library(ggplot2)
  colnames(df) <- c("group", "values")
  
  # Initialize plot
  p <- ggplot(df, aes(x = group, y = values))  
  
  # Add box
  if(isTRUE(box)){
    p <- p +
      geom_boxplot(aes(fill = group, color=group), size=1, na.rm = T, alpha = 0.4, outlier.colour = NA) + 
      geom_jitter(width = 0.2, alpha = 0.6, pch = 21, size=1, stroke=1, fill="black", color="black", na.rm = T) 
  } else {
    p <- p +
      geom_jitter(aes(fill = group, color=group), pch = 21, size = 3, width = 0.4, alpha = 0.7, na.rm = T) 
  }
  
  # Add ANOVA/Tukey results as caption
  if(add_Tukey){
   source("../functions_ANOVA_Tukey.R") 
   caption <- paste(caption, get_ANOVA_Tukey_text(df))
  }
  
  if(!is.null(pval.test)){
    # Make list of unique elements
    ele <- unique(as.character(df$group))
    # Make list of combinations (order doesn't matter) for p-values
    comb <- gtools::combinations(n = length(ele), r = 2, v = ele, repeats.allowed = F) %>% # gtools
      split(., seq(nrow(.)))
    
    # Add significance levels
    # Star height relative to bars
    vjust <- ifelse(pval.label == "p.signif", 0.5, -0.1)
    
    p <- p + stat_compare_means(
      method = pval.test, comparisons = comb, na.rm = T, vjust = vjust, hide.ns = T,
      label = pval.label, size = font_size / 2.5, bracket.size = 1
    )
  }
  # Add points, captions, etc
  p <- p +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_light() +
    labs(
      title = title,
      subtitle = subtitle,
      y = ylab,
      x = xlab,
      caption = caption
    ) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size=15),
          axis.text.x = element_text(size=25, face = "bold", colour = "black"),
          plot.caption = element_text(size=13, face = "bold", colour = "black"),
          text = element_text(size=20)
          )
  
  # Add a line for mean
  if(isTRUE(mean_line)){
    p <- p +
      stat_summary(
        fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
        position = position_dodge(width = 0.3), width = 0.25, color = "black", size = 2
      )
    
  }
  
  return(p)
}

