#' Purpose:
#' To get a stacked barplot showing the number of cells in each cc phase (y axis=count) in each celltype (x, bar)   
#' 
#' 
#' Total Cell count
# table(b2$broadcelltype)
# Mature luminal Luminal progenitor              Basal 
# 2172               3364                         1272 

# Remove all objects in environment except b2 (b2 = Seurat object)
rm(list=ls()[!ls()%in%"b2"])

# Import libraries
library(kazutils)
load_packages(c("Seurat", "ggplot2", "dplyr", "viridis"))

# Out directory
out_dir <- kazutils::create_folder("Final Plots")

# # Data 
# b2 <- readRDS("../annotated_Seurat_b2.rds")

# Colors
Phase_colors <- readRDS("../colors_Revelio_ccPhase.rds")

# New identities
Idents(b2) <- b2$broadcelltype2

# Bar plot of cell counts
celltypes <- b2$broadcelltype2 %>% unique

# Data frame of frequency of each phase in each celltype
df <- lapply(celltypes, function(ct){
  table(b2$Revelio_ccPhase[b2$broadcelltype2 == ct]) %>%
    data.frame()
})
names(df) <- celltypes
# List to one data frame
df <- do.call(rbind.data.frame, df)
# Add column for cell type
df$celltype <- kazutils::get_nth_part(rownames(df), "\\.", 1)
head(df)

# Inititate file
pdf(sprintf("%s/1_cell_number.pdf", out_dir))
# Plot
ggplot(df, aes(x = celltype, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge(), width = .8, na.rm = T) + # bars
  scale_fill_manual(name = "Phase", values = Phase_colors) +
  labs(
    title = "Number of Cells in Each Phase", # labels
    subtitle = "Phase assigned using Revelio package",
    x = "Phase",
    y = "Number of Cells"
  )  +
  geom_text(aes(label=Freq), vjust = -0.4, size=5, position = position_dodge(1))+
  theme_light() +
  lims(y=c(0,1100)) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size=15),
        axis.text.x = element_text(size=25, face = "bold", colour = "black"),
        plot.caption = element_text(size=13, face = "bold", colour = "black"),
        text = element_text(size=20)
  )
# End file
dev.off()

# Write to file
# > df
#     Var1 Freq
# LP.1 M.G1  228
# LP.2 G1.S  480
# LP.3    S  444
# LP.4   G2  933

# Counts
cnts <- table(b2$Revelio_ccPhase)
# M.G1 G1.S    S   G2 G2.M 
# 803 1111 1206 1192  495 

# Make data frame
df2 <- data.frame(CellType = c(get_nth_part(rownames(df), "\\.", 1), rep("Total", length(cnts))),
                  CellCyclePhase = c(as.character(df$Var1), names(cnts)),
                  CellCount = c(df$Freq, as.numeric(cnts)))

head(df2)

# Write data to file
write.csv(df2, file = sprintf("%s/1_cell_number_Revelio_ccPhase.csv", out_dir), row.names = F)


