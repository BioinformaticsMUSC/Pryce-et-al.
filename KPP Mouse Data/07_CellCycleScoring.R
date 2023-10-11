library(RCurl)
library(AnnotationHub)
library(ensembldb)

setwd("/Users/SuganyaSubramanian/BenKpp/")

#### Negative data #####
load("Negative/output_Relabel/09_SeuratObj_SCT_30pcs_04res_labelled.RData")

load("/Users/SuganyaSubramanian/MouseCellCycleGenes.rda")

head(g2m_genes)
head(s_genes)

seuObject_slim_nodoub_slim <- CellCycleScoring(seuObject_slim_nodoub_slim,
                                               g2m.features = g2m_genes,
                                               s.features = s_genes)

save(seuObject_slim_nodoub_slim, file ="Negative/output_Relabel/10_SeuratObj_SCT_30pcs_04res_labelled_WithCellCycle.RData")

pdf("Negative/output_Relabel/06_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuObject_slim_nodoub_slim,
        reduction = "umap",
        group.by= "Phase")
dev.off()

rm(list=ls())

##### Positive #####
load("Positive/output_Relabel/09_SeuratObj_SCT_30pcs_06res_labelled.RData")

load("/Users/SuganyaSubramanian/MouseCellCycleGenes.rda")

seuObject_slim_nodoub_slim <- CellCycleScoring(seuObject_slim_nodoub_slim,
                                               g2m.features = g2m_genes,
                                               s.features = s_genes)

save(seuObject_slim_nodoub_slim, file ="Positive/output_Relabel/10_SeuratObj_SCT_30pcs_06res_labelled_WithCellCycle.RData")

pdf("Positive/output_Relabel/06_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuObject_slim_nodoub_slim,
        reduction = "umap",
        group.by= "Phase")
dev.off()

load("Negative/output_Relabel/10_SeuratObj_SCT_30pcs_04res_labelled_WithCellCycle.RData")
negSeurat <- seuObject_slim_nodoub_slim
save(negSeurat, file ="Negative/output_Relabel/10_SeuratObj_SCT_30pcs_04res_labelled.RData_WithCellCycle.RData")
rm(seuObject_slim_nodoub_slim)

load("Positive/output_Relabel/10_SeuratObj_SCT_30pcs_06res_labelled_WithCellCycle.RData")
posSeurat <- seuObject_slim_nodoub_slim
save(posSeurat, file ="Positive/output_Relabel/10_SeuratObj_SCT_30pcs_06res_labelled_WithCellCycle.RData")
rm(seuObject_slim_nodoub_slim)
