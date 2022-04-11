#' Functions to perform ANOVA on each group
#' @param df2 Data frame with two columns in this order: 1) groups (factor/character), 2) numeric values
perform_ANOVA_Tukey <- function(df2){
  colnames(df2) <- c("group2", "vals")
  # Compute ANOVA
  df <- aov(df ~ group2, df2)
  # Extract F value
  fval <- summary(df)[[1]]$"F value"[[1]]
  # P value
  pval <- summary(df)[[1]]$"Pr(>F)"[[1]]
  # Tukey HSD
  q <- TukeyHSD(df)
  
  # Return
  c(F_value = fval, 
    p_value = pval,
    TukeyHSD = q$group2)
}


#' Get label text of values
#' @param df2 Data frame with two columns in this order: 1) groups (factor/character), 2) numeric values
get_ANOVA_Tukey_text <- function(df2, label =""){
  colnames(df2) <- c("group2", "vals")
  
  # Compute the analysis of variance
  # res.aov <- aov(df2$vals ~ df2$group2)
  res.aov <- aov(vals ~ group2, data = df2)
  # Summary of the analysis
  summary(res.aov)
  #reformat output of anova
  anova_output <- anova(res.aov)
  # get p-val
  anova_pval <-anova_output $`Pr(>F)`[1] %>% 
    formatC(format = "e", digits = 2)
  # Tukey HSD (Tukey Honest Significant Differences)
  # for performing multiple pairwise-comparison between the means of groups.
  tukey_output <- as.data.frame(TukeyHSD(res.aov)[[1]])
  tukey <- formatC(tukey_output$`p adj`, format = "e", digits = 2)
  names(tukey) <- rownames(tukey_output)
  # Return
  sprintf("%s ANOVA %s. Tukey %s", label, anova_pval, paste(names(tukey), tukey, sep = ":", collapse = ", "))
}
