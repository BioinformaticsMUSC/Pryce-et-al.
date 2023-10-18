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
source("/Users/suganyasubramanian/Shikhar/Utils.R")
setwd("/Users/SuganyaSubramanian/Guttridge_Human/")

load("output_sct_Integrated_Positive/05_SeuratObj_SCT_30pcs_05res_Slimmed.RData")

#Convert into Single Cell Experiment object
sce <- as.SingleCellExperiment(seuObject_integrated_slim)

# Doubleting by genotype
sce <- scDblFinder(sce,
                   samples="Genotype", 
                   BPPARAM=MulticoreParam(3),
                   nfeatures = 3000,
                   dims = 30,
                   dbr.sd = 1)

# add in the new column "Doublets" to make if it doublet or singlet
seuObject_integrated_slim@meta.data$Doublets <- sce$scDblFinder.class

save(seuObject_integrated_slim, file = "output_sct_Integrated_Positive/05_SeuratObj_SCT_slim_withDoublets.RData")

df <- seuObject_integrated_slim@meta.data %>% as.data.frame()
tmp <- table(df$Doublets,df$Sample) %>% as.data.frame()
tmp

pdf("output_sct_Integrated_Positive/09_GeneExp_Singlet-Doublet_Sample.pdf",width = 6,height = 6)
ggplot(tmp,aes(x=Var2, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  #scale_fill_manual(values=c("orange2", "green3")) +
  geom_text(aes(label=Freq), size=3.5) +
  ggtitle("Singlet/Doublet by Sample") +
  xlab("Sample") +
  ylab("Gene Expression") +
  theme(legend.position="none")
dev.off()

save(seuObject_integrated_slim, file = "output_sct_Integrated_Positive/06_SeuratObj_SCT_withDoublet.RData")

#Subset only singlets (remove Doublets)
seuObject_slim_nodoub <- subset(seuObject_integrated_slim, subset = Doublets == "singlet")

seuObject_slim_nodoub <- processing_seurat_sctransform(seuObject_slim_nodoub, 
                                                       vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                                       npcs = 30, 
                                                       res = c(0.3,0.4,0.5,0.6))

DefaultAssay(seuObject_slim_nodoub) <- "RNA"
seuObject_slim_nodoub <- NormalizeData(object = seuObject_slim_nodoub, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = 10000)


save(seuObject_slim_nodoub, file = "output_sct_Integrated_Positive/06_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered.RData")

seuObject_slim_nodoub_slim <- DietSeurat(seuObject_slim_nodoub, 
                                         counts = TRUE, 
                                         data = TRUE, 
                                         scale.data = FALSE,
                                         assays="RNA",
                                         dimreducs = c("pca","umap"))

save(seuObject_slim_nodoub_slim, file = "output_sct_Integrated_Positive/06_SeuratObj_SCT_30pcs_06res_NoDoublet_Reclustered_slim.RData")



pdf("output_sct_Integrated_Positive/09_UMAP_NoDoublets_Reclustered.pdf", width = 13, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Sample")
plot_grid(p1, p2)
dev.off()

library(harmony)
seuObject_slim_nodoub_withHarmony <- RunHarmony(seuObject_slim_nodoub, assay.use="RNA", group.by.vars = "Sample")
seuObject_slim_nodoub_withHarmony <- RunUMAP(seuObject_slim_nodoub_withHarmony, reduction = "harmony", dims = 1:30)
seuObject_slim_nodoub_withHarmony <- FindNeighbors(seuObject_slim_nodoub_withHarmony, reduction = "harmony", dims = 1:30) 
seuObject_slim_nodoub_withHarmony <- FindClusters(seuObject_slim_nodoub_withHarmony, resolution = c(0.3,0.4,0.5,0.6,0.7), 
                                                  algorithm = 1, n.iter = 1000,graph.name = "SCT_snn")

head(seuObject_slim_nodoub_withHarmony)

pdf("output_sct_Integrated_Positive/10_UMAP_NoDoublets_Reclustered_AfterHarmony_Res07.pdf", width = 13, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub_withHarmony, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub_withHarmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Sample")
plot_grid(p1, p2)
dev.off()

Idents(seuObject_slim_nodoub_withHarmony) <- "SCT_snn_res.0.3"
seuObject_slim_nodoub_withHarmony$seurat_clusters <- seuObject_slim_nodoub_withHarmony$SCT_snn_res.0.3

pdf("output_sct_Integrated_Positive/10_UMAP_NoDoublets_Reclustered_AfterHarmony_Res03.pdf", width = 13, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub_withHarmony, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub_withHarmony, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Sample")
plot_grid(p1, p2)
dev.off()

save(seuObject_slim_nodoub_withHarmony, file = "output_sct_Integrated_Positive/07_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered_withHarmony.RData")

#library(grDevices)
#color <- c("#17154f","#b0799a","#e69b00","#355828")
#color <- rainbow(34) # Horrible color combination

#library("randomcoloR")
#color <- distinctColorPalette(43)                    # Sample colors

color <- c("#17154f","#b0799a","#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

prop_cell <-  seuObject_slim_nodoub_withHarmony@meta.data %>%
  as.data.frame() %>%
  arrange(seurat_clusters) %>%
  group_by(seurat_clusters,Sample) %>% 
  summarise(cnt = n()) %>% 
  group_by(seurat_clusters) %>%  
  mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  group_by(seurat_clusters) %>%
  ggbarplot("seurat_clusters", "percent",
            fill = "Sample", color = "Sample", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_sct_Integrated_Positive/11_CellProportion_By_Sample.pdf", plot = prop_cell, width = 25, height = 15, units = "in", dpi = 150)

pdf("output_sct_Integrated_Positive/12_UMAP_Nodoublets_SplitBy_Sample.pdf", width = 10, height = 5)
DimPlot(seuObject_slim_nodoub_withHarmony,split.by = "Sample",ncol = 3) + theme(legend.position="none")
dev.off()
