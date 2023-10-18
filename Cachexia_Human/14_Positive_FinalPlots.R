setwd("/Users/SuganyaSubramanian/Guttridge_Human/")

load("output_Relabel_Positive/08_Seurat_Positive_Final_Relabeled.RData")

dir.create("output_Figure_Positive")
umap <- as.data.frame(Embeddings(Pos_Seurat, reduction = "umap"))
meta <- as.data.frame(Pos_Seurat@meta.data)


df <- cbind(umap,meta)%>% 
  dplyr::group_by(Cell) %>% 
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Sample) %>%
  as.data.frame()

label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  dplyr::group_by(Cell) %>% 
  dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  dplyr::left_join(label) %>%
  as.data.frame() 


#choose No. depending on the total no. of clusters
#col <- randomcoloR::distinctColorPalette(24)

#write.table(col, "output_Figure_Positive/Colors_Used.txt", sep="\t",quote=F)

#Manually edited the colors
#col <- data.table::fread("output_sct_0.6/Colors_Used_Initial_Data_edited.txt")
#col <- as.character(col$V2)
#Macro #7FE5E1
#Endo #E07647
#Epi #D2E3C3
#Erythro #DFE69C
#skeletal #7EE6B8
#pericytes #80A878
#LEC #72A1DD
#schwann #E4A6A0
#Teno #6F5CD8
#TCells #7BED43
#monocytes #DEB565
#MuSC #886E6E

#col <- c("#7FE5E1","#E07647","#D4B8E3","#DFE69C",
#        "#7EE6B8","#80A878","#72A1DD","#E4A6A0","#6F5CD8","#7BED43",
#        "#DEB565","#886E6E")


pdf("output_Figure_Positive/08_UMAP_Labelled_Alternative.pdf", width = 7, height = 6)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = Cell),size=0.2) +
  ggrepel::geom_text_repel(data = label_2, aes(label = label),
                           color = "black",
                           #fontface = 'bold',
                           segment.colour = "grey60",
                           box.padding = unit(0.25, "lines"),
                           point.padding = unit(0.5, "lines"),
                           nudge_x = .15,
                           nudge_y = 1,
                           size = 6) + 
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  theme_classic()+ #theme_void() + #this gives x axis and y axis lines
  theme(legend.position="none")
dev.off()

pdf("output_Figure_Positive/09_UMAP_Labelled_Alternative_withlegend_withaxis.pdf", width = 6, height = 5)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = Cell),size=0.1) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  #theme(legend.key.size = unit(100, 'cm'))+ to increase dot size - dint work
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic()#theme_void()
dev.off()

geno_order <- c("Pos_Control","Pos_WS_PDAC","Pos_C_PDAC")

pdf("output_Figure_Positive/09_UMAP_Labelled_Alternative_withlegend_withaxis_splitby.pdf", width = 18, height = 6)
df %>%
  mutate(Geno_ordered=fct_relevel(Sample,geno_order)) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = Cell),size=0.1) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  #theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  facet_wrap(~Geno_ordered)
dev.off()

##### Subset Macro and Recluster

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  #library(speckle)
  #library(magrittr)
  #library(broom)
  #library(muscat)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  #library(scds)
})

source("/Users/SuganyaSubramanian/Utils.R")

#Subset macrophage
pos_submacro <- subset(Pos_Seurat, subset = Cell == "Macrophage")
table(pos_submacro$Cell)

pos_submacro <- processing_seurat_sctransform(pos_submacro, 
                                             vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                             npcs = 30, 
                                             res = c(0.3))#,0.4,0.5

DefaultAssay(pos_submacro) <- "RNA"
pos_submacro <- NormalizeData(object = pos_submacro, 
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)

save(pos_submacro, file = "output_Figure_Positive/04_SubsetMacro_Reclustered.RData")

#load("output_Figure_Positive/04_SubsetMacro_Reclustered.RData")
#DimPlot(object = pos_submacro, group.by = "SCT_snn_res.0.2",reduction = "umap", 
#        label = TRUE, pt.size = 0.5) + theme(legend.position="none")

#DimPlot(object = pos_submacro, group.by = "SCT_snn_res.0.3",reduction = "umap", 
#        label = TRUE, pt.size = 0.5) + theme(legend.position="none")

#DimPlot(object = pos_submacro, group.by = "SCT_snn_res.0.4",reduction = "umap", 
#        label = TRUE, pt.size = 0.5) + theme(legend.position="none")

#Idents(pos_submacro) <- "SCT_snn_res.0.3"
#pos_submacro@meta.data$seurat_clusters <- pos_submacro@meta.data$SCT_snn_res.0.3

pdf("output_Figure_Positive/10_Macro_SubsetUMAP.pdf", width = 8, height = 4)
p1 <- DimPlot(object = pos_submacro, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = pos_submacro, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Sample")
plot_grid(p1, p2)
dev.off()

pdf("output_Figure_Positive/11_Macro_SubsetUMAP_Splitby_Genotype.pdf", width = 12, height = 4)
DimPlot(object = pos_submacro, reduction = "umap", label = TRUE,ncol = 3, pt.size = 0.5,split.by = "Sample") + theme(legend.position="none")
dev.off()

#### find markers for subMacro
all_markers_clustID <- presto::wilcoxauc(pos_submacro, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Figure_Positive/Macro_MarkersbyCluster_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

### Violin plots

#genes <- c("Ccl2","Cxcl1","Cxcl2","Cxcl5","Cxcl9","Cxcl14","Cxcl13","Ccl7","Sparc","Sparcl1")

genes <- c("MAF","MEF2C","LYVE1","TCF4","CCR2","S100A9","IL1B")

#pdf("output_Figure_Positive/12_Macro_SubsetUMAP_ViolinStacked.pdf", width = 10, height = 20)
#VlnPlot(pos_submacro, genes, stack = TRUE, sort = TRUE, flip = TRUE) +
#  theme(legend.position = "none") + ggtitle("Identity on x-axis")
#dev.off()

pdf("output_Figure_Positive/13_Macro_SubsetUMAP_Violin_1.pdf", width = 20, height = 15)
VlnPlot(pos_submacro, genes,ncol = 3,pt.size = 0)
dev.off()

### Percentage plots
#Percentages of each seurat cluster
#https://r-graph-gallery.com/128-ring-or-donut-plot.html
library(ggplot2)
library(RColorBrewer)
library(reshape2)

df <- table(Idents(pos_submacro),pos_submacro$Sample)
df <- as.data.frame(df)

df_ctrl <- df[df$Var2=="Pos_Control",]
df_ctrl$Var1 <- as.character(df_ctrl$Var1)

# Compute percentages
df_ctrl$fraction = df_ctrl$Freq / sum(df_ctrl$Freq)

# Compute the cumulative percentages (top of each rectangle)
df_ctrl$ymax = cumsum(df_ctrl$fraction)

# Compute the bottom of each rectangle
df_ctrl$ymin = c(0, head(df_ctrl$ymax, n=-1))

p1 <- ggplot(df_ctrl) +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, 
                fill=Var1), color="black") +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
  theme_void() +
  ggtitle("Control") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")# Try to remove that to see how to make a pie chart

df_cpdac <- df[df$Var2=="Pos_C_PDAC",]
df_cpdac$Var1 <- as.character(df_cpdac$Var1)

# Compute percentages
df_cpdac$fraction = df_cpdac$Freq / sum(df_cpdac$Freq)

# Compute the cumulative percentages (top of each rectangle)
df_cpdac$ymax = cumsum(df_cpdac$fraction)

# Compute the bottom of each rectangle
df_cpdac$ymin = c(0, head(df_cpdac$ymax, n=-1))

p2 <- ggplot(df_cpdac) +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, 
                fill=Var1), color="black") +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("C_PDAC") + # Try to remove that to see how to make a pie chart
  theme(plot.title = element_text(hjust = 0.5))

df_wspdac <- df[df$Var2=="Pos_WS_PDAC",]
df_wspdac$Var1 <- as.character(df_wspdac$Var1)

# Compute percentages
df_wspdac$fraction = df_wspdac$Freq / sum(df_wspdac$Freq)

# Compute the cumulative percentages (top of each rectangle)
df_wspdac$ymax = cumsum(df_wspdac$fraction)

# Compute the bottom of each rectangle
df_wspdac$ymin = c(0, head(df_wspdac$ymax, n=-1))

p3 <- ggplot(df_wspdac) +
  geom_rect(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, 
                fill=Var1), color="black") +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("WS_PDAC") + # Try to remove that to see how to make a pie chart
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")# Try to remove that to see how to make a pie chart


# Make the plot
pdf("output_Figure_Positive/14_Macro_SubsetUMAP_DonutPlot.pdf", width = 12, height = 4)
p1 + p3 + p2
dev.off()

##### final Macro UMAPS
umap <- as.data.frame(Embeddings(pos_submacro, reduction = "umap"))
meta <- as.data.frame(pos_submacro@meta.data)

df <- cbind(umap,meta)%>% 
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Sample) %>%
  as.data.frame()

label <- data.frame(seurat_clusters=levels(df$seurat_clusters),label=levels(df$seurat_clusters))

label_2 <- df %>% 
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  dplyr::left_join(label) %>%
  as.data.frame() 

pdf("output_Figure_Positive/15_Macro_FinalUMAP_withlegend_withaxis.pdf", width = 5, height = 4)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = seurat_clusters),size=0.1) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic()#theme_void()
dev.off()

#geno_order <- c("Pos_Control","Pos_WS_PDAC","Pos_C_PDAC")

pdf("output_Figure_Positive/16_Macro_FinalUMAP_withlegend_withaxis_splitby.pdf", width = 18, height = 6)
df %>%
  mutate(Geno_ordered=fct_relevel(Sample,geno_order)) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = seurat_clusters),size=0.3) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  #theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  facet_wrap(~Geno_ordered)
dev.off()























