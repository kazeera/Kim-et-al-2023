# Kim-et-al-2022
This repository contains code used to interrogate and visualize omics and other data published Kim et al. 2022, Nature Communications.

### Acronym Dictionary
Here is a list of acronyms and their respective expansions that are commonly used in the code.

**Mammary epithelial cell types/populations**
* BC = basal cells, depending on context "BC" could also refer to "Breast Cancer"
* LP = luminal progenitors
* LM, ML, LC = mature luminal cells

**Other**
* cc = cell cycle
* DDR = DNA damage response
* df = data frame or table
* DSB = double-stranded breaks
* FC = fold-change
* GDSC = external dataset, expanded is: Genomics of Drug Sensitivity in Cancer [URL](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html)
* h = human
* HK = Hyeyeon Kim, first author
* ibaq-adj lfq = iBAQ and LFQ (label-free quantitation) are measurements of protein intensity in proteomics, here it is iBAQ-adjusted by LFQ, see Methods
* m = mouse
* MEC = mammary epithelial cell(s)
* PAM50 = breast cancer subtyping method, expanded is: Prediction Analysis of Microarray 50
* PCC = Pearson correlation coefficient
* ssgsea = Computational method from GSVA package, single-sample gene set enrichment analysis



### R Packages Used
R Package | Version | Use | Reference
:- |:- |:- |:- |
R | 4.1.0 | Computing language | R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. 
Seurat | 4.0.2 | Single-cell (sc) object | Hao and Hao et al (2020). Integrated analysis of multimodal single-cell data. bioRxiv. [URL](https://satijalab.org/seurat/)
Revelio | 0.1.0 | Sc cell cycle annotation | Schwabe, Daniel, et al (2020). The transcriptome dynamics of single cells during the cell cycle. Molecular systems biology 16.11. [URL](https://www.embopress.org/doi/full/10.15252/msb.20209946)
plyr | 1.8.6 | Mapvalues | Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. [URL](http://www.jstatsoft.org/v40/i01/)
viridis | 0.6.1 | Sc color palette | Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2021). Rvision - Colorblind-Friendly Color Maps for R. R package version 0.6.1.
HGNChelper | 0.8.1 | Check symbols for Revelio cell cycle markers | Levi Waldron and Markus Riester (2019). HGNChelper: Identify and Correct Invalid HGNC Human Gene Symbols and MGI Mouse Gene Symbols. R package version 0.8.1. [URL](https://f1000research.com/articles/9-1493)
ggplot2 | 3.3.5 | All plots | H. Wickham (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. [URL](https://www.springer.com/gp/book/9780387981413)
dplyr | 1.0.6 | Data manipulation | Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.6. 
escape | 1.2.0 | Sc gene set enrichment | Nick Borcherding and Jared Andrews (2021). escape: Easy single cell analysis platform for enrichment. R package version 1.2.0.
GSEABase | 1.54.0 | Sc GSEA data structure | Martin Morgan, Seth Falcon and Robert Gentleman (2021). GSEABase: Gene set enrichment data structures and methods. R package version 1.54.0.
ggridges | 0.5.3 | Ridge plot | Claus O. Wilke (2021). ggridges: Ridgeline Plots in 'ggplot2'. R package version 0.5.3.
GSVA | 1.40.0 | ssGSEA | Hänzelmann, S., Castelo, R. and Guinney, A (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics, 14:7. [URL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)
ComplexHeatmap | 2.8.0 | PDX heatmap | Gu, Z (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. [URL](https://academic.oup.com/bioinformatics/article/32/18/2847/1743594)
circlize | 0.4.12 | Used to make ComplexHeatmap annotation | Gu, Z (2014). circlize implements and enhances circular visualization in R. Bioinformatics. [URL](https://academic.oup.com/bioinformatics/article/30/19/2811/2422259)
ggpubr | 0.4.0 | Stats on ggplot | Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0.
reshape2 | 1.4.4 | Reformat dfs | Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. [URL](http://www.jstatsoft.org/v21/i12/)
pheatmap | 1.0.12 | All other heatmaps | Raivo Kolde (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
RColorBrewer | 1.1.2 | RdBu palette | Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1.2.
kazutils | | Utility functions | Kazeera Aliar (2021). [URL](https://github.com/kazeera/kazutils)
