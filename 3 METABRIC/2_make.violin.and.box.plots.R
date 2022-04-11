#' Purpose:
#' Make box plots comparing signature enrichments in PAM50 subtypes and perform significance testing

# Save plots
library("optparse")
library(dplyr)
library(foreach)
library(ggplot2)
library(RColorBrewer)

# Load data 
# ssGSEAs.files <- list.files(pattern=".csv")
signatures <- colnames(ssGSEA.subtype.list[[1]][2:4])

#Reading ssGSEAs files Function
ReadssGSEAFunc <- function(ssGSEAs.files, subtype) {
  ssGSEA = read.csv(file.path(subtypes.path,ssGSEAs.files)[subtype])
}

# signatures and subtypes data frame
ssGSEA.subtype = do.call(rbind.data.frame, ssGSEA.subtype.list)

ssGSEA.subtype = ssGSEA.subtype %>% select(
  subtype = Pam50...Claudin.low.subtype,
  Basal = Basal,
  `Luminal Mature` = `Luminal Mature`,
  `Luminal Progenitor` = `Luminal Progenitor`
)
# Stats functions ====================================
# ttest function
ttestFunc <- function (ssGSEA.subtype, signature) {
  signature.values = ssGSEA.subtype[,signature]
  subtype = ssGSEA.subtype$subtype
  
  test = pairwise.t.test(
    signature.values, subtype, 
    p.adjust.method = "fdr"
  )
  print(signature)
  print(test$method)
  print(test$p.value)
}

# ANOVA Func
anovaTestFunc <- function(ssGSEA.subtype, signature){
  # run 1 way ANOVA across all values
  res.aov <- aov(formula = ssGSEA.subtype[,signature] ~ ssGSEA.subtype$subtype)
  
  return(res.aov)
}
# Plotting functions =================================
# Violin Plot Function
drawViolinPlot <- function(ssGSEA.subtype, signature, aovpvalue){
  # Subset data based on signature of interest 
  dfAux = data.frame(
    subtype = ssGSEA.subtype$subtype,
    signature = ssGSEA.subtype[,signature]
  )  
  # Sort x axis subtypes so it's consistent (based on defined order)
  dfAux$subtype <- factor( dfAux$subtype,levels = 
                             c("claudin-low", "Basal", "LumA", "LumB", "Her2", "NC", "Normal"),ordered = TRUE)
  # Plot 
  g <- ggplot(dfAux, aes(x = subtype, y = signature,)) + 
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill=NA)+#, outlier.shape = NA) +
    # stat_summary(fun.y=median, geom="line", size=2, color="red")+
    geom_jitter(width = 0.2)+
    # geom_point(pch = 21, position = position_jitterdodge())+
    # scale_y_continuous(name = "ssGSEA Enrichment Score", breaks = seq(-2, 2, .2))+
    ylim(-0.42, 0.69)+
    # scale_fill_manual(values = unname(c(color = brewer.pal(length(levels(dfAux$subtype)), "Greys"))))+
    xlab("Breast Cancer Subtype")+
    ylab("ssGSEA Enrichment Score")+
    # geom_point()+
    theme_classic() 
  return(g)
}

# Boxplot Function
drawBoxPlot <- function(ssGSEA.subtype, signature, aovpvalue){
  dfAux = data.frame(
    subtype = ssGSEA.subtype$subtype,
    signature = ssGSEA.subtype[,signature]
  )  
  
  # x <- unlist(lapply(split(dfAux$signature,dfAux$subtype),median))
  # y <- sort(x, decreasing = TRUE, na.last = TRUE)
  # 
  # dfAux$subtype <- factor( dfAux$subtype,levels = names(y),ordered = TRUE)
  dfAux$subtype <- factor( dfAux$subtype,levels = 
                             c("claudin-low", "Basal", "LumA", "LumB", "Her2", "NC", "Normal"),ordered = TRUE)
  
  g <-ggplot(dfAux, aes(x = subtype, y = signature)) + 
    geom_boxplot(data = dfAux ) +
    # geom_jitter(width=0.1)+
    # scale_y_continuous(name = "ssGSEA Enrichment Score", breaks = seq(-2, 2, .2))+
    ylim(-0.42, 0.69)+
    # ylab("ssGSEA Enrichment Score")+
    xlab("Breast Cancer Subtype")+
    scale_x_discrete(labels=c("claudin-low"="CL", "Basal"="BC", "LumA"="LA", "LumB"="LB", "Her2"="H2", "NC"="NC", "Normal"="No"))+
    # geom_point()+
    theme_classic() +
    ggtitle(signature, subtitle = sprintf("ANOVA : %e", aovpvalue)) +
    ylab("ssGSEA Enrichment Score")
  return(g)
}

# Main loop to make plots ==============================
for (signature in signatures){
  # get anova results from ssGSEA across signatures
  res.aov <- anovaTestFunc(ssGSEA.subtype, signature)
  #reformat output of anova
  anova_output <- anova(res.aov)
  # get p-val
  anova_pval <-anova_output $`Pr(>F)`[1]
  
  # Make violin plot
  violinPlot <- drawViolinPlot(ssGSEA.subtype, signature, anova_pval)
  ggsave(sprintf("%s/violinplot_BC.subtype_human.sig.FC%s_%s.png", output_folder, FC, signature), violinPlot)
  ggsave(sprintf("%s/violinplot_BC.subtype_human.sig.FC%s_%s.svg", output_folder, FC, signature), violinPlot)
  ggsave(sprintf("%s/violinplot_BC.subtype_human.sig.FC%s_%s.pdf", output_folder, FC, signature), violinPlot)
  
  # Make box plot
  boxPlot <- drawBoxPlot(ssGSEA.subtype, signature, anova_pval)
  ggsave(sprintf("%s/boxplot_BC.subtype_human.sig.FC%s_%s.png", output_folder, FC, signature), boxPlot)
  ggsave(sprintf("%s/boxplot_BC.subtype_human.sig.FC%s_%s.pdf", output_folder, FC, signature), boxPlot)
  
  #reformat output of anova
  anova_output <- anova(res.aov)
  # get p-val
  anova_pval <-anova_output $`Pr(>F)`[1]
  # run TukeyHSD test
  tukey_output <- as.data.frame(TukeyHSD(res.aov)[[1]])
  
  # Write stats to file
  sink(sprintf("%s/anova.tukey_human.sig.FC%s_%s.txt", output_folder, FC, signature))
  print(res.aov)
  print(anova_output)
  print(tukey_output)
  
  sink()
  x<-ttestFunc(ssGSEA.subtype, signature)
  sink(sprintf("%s/t.tests_human.sig.FC%s_%s.txt", output_folder, FC, signature))
  print(x)
  sink()
}
