setwd("/Users/SuganyaSubramanian/BenKpp/")
load("output_Relabel_Negative/11_FinalRelabel.RData")

umap <- as.data.frame(Embeddings(negSeurat, reduction = "umap"))
meta <- as.data.frame(negSeurat@meta.data)

df <- cbind(umap,meta)%>% 
  dplyr::group_by(Cell) %>% 
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Genotype) %>%
  as.data.frame()

label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  dplyr::group_by(Cell) %>% 
  dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  dplyr::left_join(label) %>%
  as.data.frame() 

#choose No. depending on the total no. of clusters

pdf("output_Relabel_Negative/08_UMAP_Labelled_Alternative.pdf", width = 7, height = 6)
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

pdf("output_Relabel_Negative/09_UMAP_Labelled_Alternative_withlegend_withaxis.pdf", width = 6, height = 5)
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

pdf("output_Relabel_Negative/09_UMAP_Labelled_Alternative_withlegend_withaxis_splitby.pdf", width = 12, height = 6)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = Cell),size=0.1) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  #theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  facet_wrap(~Genotype)
dev.off()

##### Subset Faps and Recluster

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

#Subset Faps
neg_subfaps <- subset(negSeurat, subset = Cell == "Faps")
table(neg_subfaps$Cell)

neg_subfaps <- processing_seurat_sctransform(neg_subfaps, 
                                                       vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                                       npcs = 30, 
                                                       res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

DefaultAssay(neg_subfaps) <- "RNA"
neg_subfaps <- NormalizeData(object = neg_subfaps, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = 10000)

save(neg_subfaps, file = "output_Relabel_Negative/04_SubsetFaps_Reclustered.RData")
load("output_Relabel_Negative/04_SubsetFaps_Reclustered.RData")
DimPlot(object = neg_subfaps, group.by = "SCT_snn_res.0.2",reduction = "umap", 
  label = TRUE, pt.size = 0.5) + theme(legend.position="none")

DimPlot(object = neg_subfaps, group.by = "SCT_snn_res.0.3",reduction = "umap", 
  label = TRUE, pt.size = 0.5) + theme(legend.position="none")

DimPlot(object = neg_subfaps, group.by = "SCT_snn_res.0.4",reduction = "umap", 
  label = TRUE, pt.size = 0.5) + theme(legend.position="none")

Idents(neg_subfaps) <- "SCT_snn_res.0.4"
neg_subfaps@meta.data$seurat_clusters <- neg_subfaps@meta.data$SCT_snn_res.0.4

pdf("output_Relabel_Negative/10_Faps_SubsetUMAP.pdf", width = 10, height = 6)
p1 <- DimPlot(object = neg_subfaps, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = neg_subfaps, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel_Negative/11_Faps_SubsetUMAP_Splitby_Genotype.pdf", width = 12, height = 7)
DimPlot(object = neg_subfaps, reduction = "umap", label = TRUE,ncol = 2, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

#### find markers for subfaps
all_markers_clustID <- presto::wilcoxauc(neg_subfaps, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Relabel_Negative/Faps_MarkersbyCluster_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

### Violin plots

genes <- c("Ccl2","Cxcl1","Cxcl2","Cxcl5","Cxcl9","Cxcl14","Cxcl13","Ccl7","Sparc","Sparcl1")

pdf("output_Relabel_Negative/12_Faps_SubsetUMAP_ViolinStacked.pdf", width = 10, height = 20)
VlnPlot(neg_subfaps, genes, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()

pdf("output_Relabel_Negative/13_Faps_SubsetUMAP_Violin_1.pdf", width = 20, height = 15)
VlnPlot(neg_subfaps, genes,ncol = 3,pt.size = 0)
dev.off()

#requested on Oct10
genes <- c("Mmp3", "Plau", "Tgif1", "Myc")
pdf("output_Relabel_Negative/13_Faps_SubsetUMAP_Violin_2.pdf", width = 20, height = 5)
VlnPlot(neg_subfaps, genes,ncol = 4,pt.size = 0)
dev.off()

### Percentage plots
#Percentages of each seurat cluster
#https://r-graph-gallery.com/128-ring-or-donut-plot.html
library(ggplot2)
library(RColorBrewer)
library(reshape2)

df <- table(Idents(neg_subfaps),neg_subfaps$Genotype)
df <- as.data.frame(df)

df_ctrl <- df[df$Var2=="CtrlNeg",]
df_ctrl$Var1 <- as.character(df_ctrl$Var1)

# Compute percentages
df_ctrl$fraction = df_ctrl$Freq / sum(df_ctrl$Freq)

# Compute the cumulative percentages (top of each rectangle)
df_ctrl$ymax = cumsum(df_ctrl$fraction)

# Compute the bottom of each rectangle
df_ctrl$ymin = c(0, head(df_ctrl$ymax, n=-1))

p1 <- ggplot(df_ctrl, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
     geom_rect() +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_void() +
     ggtitle("Control") +
     theme(plot.title = element_text(hjust = 0.5)) +
     theme(legend.position = "none")# Try to remove that to see how to make a pie chart

df_kpp <- df[df$Var2=="KppNeg",]
df_kpp$Var1 <- as.character(df_kpp$Var1)

# Compute percentages
df_kpp$fraction = df_kpp$Freq / sum(df_kpp$Freq)

# Compute the cumulative percentages (top of each rectangle)
df_kpp$ymax = cumsum(df_kpp$fraction)

# Compute the bottom of each rectangle
df_kpp$ymin = c(0, head(df_kpp$ymax, n=-1))

p2 <- ggplot(df_kpp, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
     geom_rect() +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) +
     theme_void() +
     ggtitle("Kpp") + # Try to remove that to see how to make a pie chart
     theme(plot.title = element_text(hjust = 0.5))

# Make the plot
pdf("output_Relabel_Negative/14_Faps_SubsetUMAP_DonutPlot.pdf", width = 8, height = 4)
p1 + p2
dev.off()

##### final FAPS UMAPS
umap <- as.data.frame(Embeddings(neg_subfaps, reduction = "umap"))
meta <- as.data.frame(neg_subfaps@meta.data)

df <- cbind(umap,meta)%>% 
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Genotype) %>%
  as.data.frame()

label <- data.frame(seurat_clusters=levels(df$seurat_clusters),label=levels(df$seurat_clusters))

label_2 <- df %>% 
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  dplyr::left_join(label) %>%
  as.data.frame() 

pdf("output_Relabel_Negative/15_FAPs_FinalUMAP_withlegend_withaxis.pdf", width = 5, height = 4)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = seurat_clusters),size=0.1) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic()#theme_void()
dev.off()

pdf("output_Relabel_Negative/15_FAPs_FinalUMAP_withlegend_withaxis_splitby.pdf", width = 10, height = 4)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = seurat_clusters),size=0.3) +
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  #scale_colour_manual(values = col)+
  #theme_void()+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  facet_wrap(~Genotype)
dev.off()






















