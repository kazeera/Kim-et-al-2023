#' Purpose:
#' To perform and visualize linear regressions of ssGSEA score vs BMN angle.
#' 
#' Notes:
#' - point-shapes: claudin-low AKA normal (X), basal-like, other
#' - label R and R^2 
#' 
#' Credit: Cescon group for data, Arvind Mer for computing angle

# Import R packages
library(kazutils) # devtools::install_github("kazeera/kazutils")
load_packages(c("openxlsx", "dplyr", "kazcolors", "ggpubr","VennDiagram", "reshape2", "ggplot2"))
source("functions_reg_plot.R")

# Output folder
out_dir <- create_folder("Plots")

# Read in data file
df <- read.csv("merged_ann_ssgsea_table.csv",stringsAsFactors = F)

# Reformat
rownames(df) <- make.unique(as.character(df$Patient.ID))
cts <- colnames(df)[grep("ssGSEA$" , colnames(df))]

## REGRESSION
# Colors
cell_type_colors <-  c(BC = "red", LM="darkblue", LP="deepskyblue3")
PAM50_colors <- c(LumA="mediumpurple1", LumB="mediumpurple4", Basal="seagreen1", Normal="seagreen", HER2="goldenrod1", Unknown="gray")

# Point shapes
df$PAM50_2 <- df$PAM50 
df$PAM50_2[! df$PAM50_2 %in% c("Basal", "Normal")] <- "Other"

# Get max and min for axes so we can post
Y_MAX <- ceiling(max(df$angle)) +20
Y_MIN <- floor(min(df$angle))
PT_SIZE <- 2.5
FONT_SIZE <- 15

# Comparisons (input columns of df)
comps <- c("BC_ssGSEA", "LP_ssGSEA", "LM_ssGSEA")

## Individual regressions with equation
pdf(sprintf("%s/2_%s_regression_ssGSEA.pdf", out_dir, format(Sys.Date(), "%Y%m%d")))
corr_method = "pearson"

for(column in comps){
        # Get column values
        celltype <- get_nth_part(column, "_", 1)
        df$value <- df[,column]
        # x position of equation on plot
        X_eqn <- (min(df$value)+max(df$value))/2
        # Plot
        p <- ggplot(data = df, aes(y=angle,x=value)) +
                geom_smooth(method="lm", size=1,se = T, fill = cell_type_colors[celltype], color="black", alpha=0.5) +
                geom_point(size = PT_SIZE, stroke = 1, alpha=0.9, color="black", aes(color=PAM50, shape=PAM50_2))+
                scale_shape_manual(values=c(Basal=2,Other=1, Normal=4))+
                # scale_fill_manual(values = cell_color)+ 
                scale_color_manual(name="PAM50", values=PAM50_colors)+
                # scale_shape_manual(values = hormone_status_shapes)+
                labs(title = sprintf("%s", celltype),
                     subtitle = "Linear regression of BMN angle vs Enrichment Score",
                     y = "Angle",
                     x = column
                )+ scale_y_continuous(limits=c(Y_MIN,Y_MAX), labels=scaleFUN)+
                geom_text(x = X_eqn, y = Y_MAX, label = lm_eqn(df), size = 3.5, parse = TRUE)+ #stats of lm+
                theme_classic()+
                theme(axis.title = element_text(size=FONT_SIZE),
                      axis.text = element_text(color="black", size=FONT_SIZE),
                      axis.ticks.length = unit(0.25, "cm"), 
                      legend.position = "bottom",
                      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), # increase space between x axis title and labels
                      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
        print(p)      
}
dev.off()
# Save workspace
save.image("workspace.Rdata")

