suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  #library(scds)
})

setwd("/Users/SuganyaSubramanian/BenKpp/")

load("output_sct_Integrated/05_SeuratObj_SCT_30pcs_05res_Slimmed.RData")

source("/Users/SuganyaSubramanian/Utils.R")

#Convert into Single Cell Experiment object
sce <- as.SingleCellExperiment(seuObject_integrated_slim)

# Doubleting by genotype
sce <- scDblFinder(sce,
                   samples="Genotype", 
                   #BPPARAM=MulticoreParam(3),
                   nfeatures = 3000,
                   dims = 30,
                   dbr.sd = 1)

# add in the new column "Doublets" to make if it doublet or singlet
seuObject_integrated_slim@meta.data$Doublets <- sce$scDblFinder.class

save(seuObject_integrated_slim, file = "output_sct_Integrated/05_SeuratObj_SCT_withDoublets.RData")

df <- seuObject_integrated_slim@meta.data %>% as.data.frame()
tmp <- table(df$Doublets,df$Genotype) %>% as.data.frame()

pdf("output_sct_Integrated/09_GeneExp_Singlet-Doublet.pdf")
ggplot(tmp,aes(x=Var2, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("orange2", "green3")) +
  geom_text(aes(label=Freq), vjust=3.0,size=3.5) +
  ggtitle("Singlet/Doublet by Genotype") +
  xlab("Sample") +
  ylab("Gene Expression") +
  theme(legend.position="none")
dev.off()

#Subset only singlets (remove Doublets)
seuObject_slim_nodoub <- subset(seuObject_integrated_slim, subset = Doublets == "singlet")

seuObject_slim_nodoub <- processing_seurat_sctransform(seuObject_slim_nodoub, 
                                                       vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                                       npcs = 30, 
                                                       res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3))

DefaultAssay(seuObject_slim_nodoub) <- "RNA"
seuObject_slim_nodoub <- NormalizeData(object = seuObject_slim_nodoub, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = 10000)

save(seuObject_slim_nodoub, file = "output_sct_Integrated/06_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered.RData")

seuObject_slim_nodoub_slim <- DietSeurat(seuObject_slim_nodoub, 
                                         counts = TRUE, 
                                         data = TRUE, 
                                         scale.data = FALSE,
                                         assays="RNA",
                                         dimreducs = c("pca","umap"))
#choose Resolution 
#0.4

Idents(seuObject_slim_nodoub_slim) <- "SCT_snn_res.0.4"
seuObject_slim_nodoub_slim@meta.data$seurat_clusters <- seuObject_slim_nodoub_slim@meta.data$SCT_snn_res.0.4

pdf("output_sct_Integrated/10_UMAP_NoDoublet.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()


save(seuObject_slim_nodoub_slim, file = "output_sct_Integrated/07_SeuObject_SCT_30pcs_allres_NoDoublet_slim.RData")

load("/Users/SuganyaSubramanian/BenKpp/output_sct_Integrated/07_SeuObject_SCT_30pcs_allres_NoDoublet_slim.RData")
seuObject_slim_nodoub_slim

all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub_slim, 'SCT_snn_res.0.4', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_sct_Integrated/Res04_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")
