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

setwd("/Users/SuganyaSubramanian/BenKpp/")

load("Negative/output_Relabel/10_SeuratObj_SCT_30pcs_04res_labelled_WithCellCycle.RData")

table(negSeurat$Cell)

###2. Load Labels data
Labels <- read.table("Negative/output_sct_Integrated/New_Labels.txt",header=T,sep="\t")


#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

current.cluster.ids
new.cluster.ids

table(Idents(negSeurat))

##Step 3: Rename cluster with new Labels
negSeurat@active.ident <- plyr::mapvalues(x = negSeurat@active.ident, 
                                                           from = current.cluster.ids, 
                                                           to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
negSeurat@meta.data$Cell <- negSeurat@active.ident

table(negSeurat@meta.data$Cell)
Faps  Endothelial   Epithelial Erythrocytes     Skeletal    Pericytes          LEC 
       12451         7289         3419         4577         3960         1573          743 
     Schwann    Tenocytes      T_Cells    Monocytes         MuSC 
         477          376          249          140           98 

##Save new Seurat data
save(negSeurat,file = "Negative/output_Relabel/11_FinalRelabel.RData")

pdf("Negative/output_Relabel/01_UMAP_Labelled.pdf", width = 6, height = 7)
DimPlot(object = negSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("Negative/output_Relabel/02_UMAP_Labelled_Groupby.pdf", width = 12, height = 7)
p1 <- DimPlot(object = negSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = negSeurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Negative/output_Relabel/03_UMAP_Labelled_withclusters.pdf", width = 12, height = 7)
Idents(negSeurat) <- "seurat_clusters"
p1 <- DimPlot(object = negSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(negSeurat) <- "Cell"
p2 <- DimPlot(object = negSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("Negative/output_Relabel/04_UMAP_Labelled_Splitby_Genotype.pdf", width = 12, height = 7)
DimPlot(object = negSeurat, reduction = "umap", label = TRUE,ncol = 2, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Faps","Skeletal","MuSC","Endothelial","Epithelial",
	"Pericytes","Schwann","Tenocytes","LEC", "T_Cells","Monocytes","Erythrocytes")
Idents(negSeurat) <- factor(Idents(negSeurat), levels = order) # reorder the factor based on cell type

markers <- c("Gsn", #"Col1a1",#Faps
             "Tnnc2","Pax7",#skeletal, MuSC
             "Pecam1","Lgals7",#"Lor", #endo, Epi"Epcam",
             "Acta2",#pericytes
             "Plp1",#Schwann
             "Thbs4", #Tenocytes "Mpz", #Tenocytes
             "Ccl21a", #LEC
             "Cd3d",#tcell
             "Fcer1g", # monocytes
             "Hbb-bs") #Tumor("Rps15a",), Erythrocytes


pdf("Negative/output_Relabel/07_DotPlot_Celltype_Markers.pdf", width = 10, height = 6)
DotPlot(negSeurat, features = markers) + rotate_x_text(angle = 45)
dev.off()