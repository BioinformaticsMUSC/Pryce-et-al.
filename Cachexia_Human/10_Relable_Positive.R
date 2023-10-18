suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(speckle)
  library(magrittr)
  library(broom)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(RColorBrewer)
  library(MetBrewer)
  library(scales)
  library(ComplexHeatmap)
  library(randomcoloR)
  library(png)#to include pic in ggplot
  library(patchwork)
})

setwd("/Users/SuganyaSubramanian/Guttridge_Human/")
#Step:1 load Seurat Object
load("output_sct_Integrated_Positive/08_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered_withHarmony_withCellScoring.RData")

#Select needed colummns only
Pos_Seurat <- seuObject_slim_nodoub_withHarmony

Pos_Seurat@meta.data <- Pos_Seurat@meta.data %>% 
  dplyr::select(orig.ident, nUMI, nGene, Sample, pMito, pRibo, Genotype,
                nCount_RNA, nFeature_RNA, nCount_SCT, nFeature_SCT,
                seurat_clusters,SCT_snn_res.0.3,Phase)

dir.create("output_Relabel_Positive")

###2. Load Labels data
Labels <- read.table("output_sct_Integrated_Positive/Labels_Cluster.txt",header=T,sep="\t")

#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

current.cluster.ids
new.cluster.ids

table(Pos_Seurat@active.ident)

##Step 3: Rename cluster with new Labels
Pos_Seurat@active.ident <- plyr::mapvalues(x = Pos_Seurat@active.ident, 
                                           from = current.cluster.ids, 
                                           to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
Pos_Seurat@meta.data$Cell <- Pos_Seurat@active.ident

##Save new Seurat data
save(Pos_Seurat,file = "output_Relabel_Positive/08_Seurat_Positive_Final_Relabeled.RData")

pdf("output_Relabel_Positive/01_UMAP_Labelled.pdf", width = 6, height = 5)
DimPlot(object = Pos_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("output_Relabel_Positive/02_UMAP_Labelled_Groupby.pdf", width = 12, height = 5)
p1 <- DimPlot(object = Pos_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Pos_Seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Sample")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel_Positive/03_UMAP_Labelled_withclusters.pdf", width = 12, height = 5)
Idents(Pos_Seurat) <- "seurat_clusters"
p1 <- DimPlot(object = Pos_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(Pos_Seurat) <- "Cell"
p2 <- DimPlot(object = Pos_Seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel_Positive/04_UMAP_Labelled_Splitby_Genotype.pdf", width = 15, height = 5)
DimPlot(object = Pos_Seurat, reduction = "umap", label = TRUE,ncol = 3, pt.size = 0.5,split.by = "Sample") + theme(legend.position="none")
dev.off()


### 2 Cell proportion ####

#color <- c("#17154f","#b0799a")
color <- c("#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

prop_cell <-  Pos_Seurat@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Cell) %>%
  dplyr::group_by(Cell,Sample) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Cell) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Cell) %>%
  ggbarplot("Cell", "percent",
            fill = "Sample", color = "Sample", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="right") +
  xlab("")

ggsave("output_Relabel_Positive/06_CellProportion_By_Sample.pdf", plot = prop_cell, width = 12, height = 6, units = "in", dpi = 150)


# 03 Cell Count
cnt <- table(Pos_Seurat$Cell, Pos_Seurat$Sample)
openxlsx::write.xlsx(cnt,file = "output_Relabel_Positive/Cell_Count.xlsx", colNames = TRUE, borders = "columns")

#04 Dot plot
order <- c("FAPs","Skeletal","Pericytes","Endothelial",
           "T_Cells","Macrophage","Neutrophils","Mastcells","Erythroblast")

Idents(Pos_Seurat) <- factor(Idents(Pos_Seurat), levels = order) # reorder the factor based on cell type

markers <- c("PDGFRA","ACTA1","ACTA2","PECAM1","CD3G","CD14","S100A9","TPSAB1","HBB")#"APOC1",

pdf("output_Relabel_Positive/07_DotPlot_Celltype_Markers.pdf", width = 10, height = 6)
DotPlot(Pos_Seurat, features = markers) + rotate_x_text(angle = 45)
dev.off()

#Presto Markers by Cell type
all_markers_clustID <- presto::wilcoxauc(Pos_Seurat, 'Cell', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers_clustID.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05)

openxlsx::write.xlsx(all_markers_clustID.Sign, file = "output_Relabel_Positive/Presto_Markers_filtered_padjlt05_CellType.xlsx", colNames = TRUE, borders = "columns")



