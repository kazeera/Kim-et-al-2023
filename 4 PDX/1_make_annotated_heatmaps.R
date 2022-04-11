#' Purpose:
#' To create a heatmap for PDX RNA expression with annotations for BMN and PAM50 subtype
#' 
#' Credit: Cescon group for data, Arvind Mer for computing angle

# Import R packages
# devtools::install_github("kazeera/kazutils")
library(kazutils)
load_packages(c("openxlsx", "dplyr", "kazcolors", "circlize", "ComplexHeatmap", "pheatmap", "RColorBrewer"))
source("make_heatmap.R") 

# Read in data file
df <- read.csv("merged_ann_ssgsea_table.csv", stringsAsFactors = F)
rownames(df) <- make.unique(df$Patient.ID)

# Out folder
output_folder <- kazutils::create_folder("Plots")

# Colors to annotate BC subtypes and "class" of cell types
PAM50_colors <- c(LumA="mediumpurple1", LumB="mediumpurple4", Basal="seagreen1", Normal="seagreen", HER2="goldenrod1", Unknown="gray")

# Make a heatmap with an annotation for BMN response   
make_heatmap(df[,grep("ssGSEA$", colnames(df))], 
             label2=paste(nrow(df), "PDX_BMN_response", sep="_"), 
             ann_row = df[,c("Main.Response", "PAM50")]
)

# For each PDX, get  a table that looks like this:
#      .      Freq
#       CR      1
#       PD      1
freq <- df$Details %>%
        gsub(";|/", ",", .) %>%
        gsub(" ", "", .) %>%
        gsub("one|mostly|some", "", .) %>%
        lapply(function(x) strsplit(x, ",") %>% unlist %>% table %>% data.frame)

# Rename column to show PDX id
# Response REF023
#       CR      1
#       PD      1
for (i in  1:length(freq)){
        colnames(freq[[i]]) <- c("Response", rownames(df)[i])
        
}

# List to data frame, where each row is a response (CR, PR, etc) and under each PDX is a sum of replicates
freq2 = bind_rows(freq) %>%
        group_by(Response) %>%
        summarise_each(funs(sum(., na.rm = T))) %>%
        # t %>%
        data.frame()

unique(freq2$Response)
# colnames(freq2) <- freq2[1,]
# freq2 <- freq2[-1,]
rownames(freq2) <- freq2$Response
freq2$Response <- NULL

# Define colors for response (same order as freq2)
colors_response <- c(CR="darkgreen", PD="pink", SD="magenta", PR="green")

# Make matrix
mat <- as.matrix(df[,c("BC_ssGSEA", "LP_ssGSEA", "LM_ssGSEA")])

# BMN angle color gradient
col_fun = circlize::colorRamp2(c(-10, 50, 164), c("orangered4", "white", "darkgreen"))
# col_fun = colorRamp2(c(min(df$angle), mean(df$angle), max(df$angle)), c("orangered4", "white", "darkgreen"))
# brew_pal <- "PuRd"
# max_n <- brewer.pal.info$maxcolors[grep(brew_pal, rownames(brewer.pal.info))]
# angle_colors <- colorRampPalette(brewer.pal(n = max_n, name = brew_pal))(50)

# Make annotations
row_ha2 <- rowAnnotation(BMN_angle = df$angle,
                         PAM50 = df$PAM50,
                         col = list(PAM50 = PAM50_colors, BMN_angle = col_fun, BMN_Response = colors_response))
# dev.off()
# Unsupervised hierarchal clustering
# win.graph()
row_ha <- rowAnnotation(Response = row_anno_barplot(t(freq2), gp=gpar(fill=colors_response)),
                        # Main = df$Main.Response,
                        col = list(Response = colors_response))

# Save to file
pdf(sprintf("1_%s_Cescon PDX heatmap.pdf", format(Sys.Date(), "%Y%m%d")))
# win.graph()
Heatmap(mat, 
        name = "ssGSEA score",
        cluster_rows = T,
        show_row_names = T,
        heatmap_width = unit(10, "cm"), 
        right_annotation = row_ha,
        left_annotation = row_ha2,
        row_names_gp = gpar(col = "black", fontsize = 7))

dev.off()

mat_s <- apply(mat, 2, function(x) scales::rescale(x, c(-2,2)))

# Save to file
pdf(sprintf("1_%s_Cescon PDX heatmap_scaled.pdf", format(Sys.Date(), "%Y%m%d")))
# win.graph()
Heatmap(mat_s, 
        name = "ssGSEA score",
        cluster_rows = T,
        show_row_names = T,
        heatmap_width = unit(10, "cm"), 
        right_annotation = row_ha,
        left_annotation = row_ha2,
        row_names_gp = gpar(col = "black", fontsize = 7))

dev.off()

win.graph()
scales::show_col(colors_response[c("CR", "PR", "SD", "PD")], ncol = 1, labels = F)# names(colors_response)) 
## Todo colors for stacked bar plot no leegend ?? 
colors_response
# See that the colors match up
rownames(freq2)
# [1] "CR" "PD" "SD" "PR"
freq2[,"REF016"]
# [1] 1 1 1 0
freq2[,"NOTCH06"]
# [1] 0 1 1 0

names(freq2) <- gsub("X", "", names(freq2))
saveRDS(freq2, "1_freq_BMN_qualitative_response.rds")
