library(RCurl)
library(AnnotationHub)
library(ensembldb)

setwd("/Users/SuganyaSubramanian/BenKpp/")

#### Negative data #####
load("output_sct_Integrated_Negative/07_SeuObject_SCT_30pcs_04res_NoDoublet_slim.RData")

load("/Users/SuganyaSubramanian/MouseCellCycleGenes.rda")

head(g2m_genes)
head(s_genes)

seuObject_slim_nodoub_slim <- CellCycleScoring(seuObject_slim_nodoub_slim,
                                               g2m.features = g2m_genes,
                                               s.features = s_genes)

save(seuObject_slim_nodoub_slim, file ="output_sct_Integrated_Negative/10_SeuratObj_SCT_30pcs_04res_labelled_WithCellCycle.RData")

pdf("output_sct_Integrated_Negative/06_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuObject_slim_nodoub_slim,
        reduction = "umap",
        group.by= "Phase")
dev.off()

negSeurat <- seuObject_slim_nodoub_slim
save(negSeurat, file ="output_sct_Integrated_Negative/10_SeuratObj_SCT_30pcs_04res_labelled.RData_WithCellCycle.RData")
rm(seuObject_slim_nodoub_slim)

##### Positive #####
load("output_sct_Integrated_Positive/07_SeuObject_SCT_30pcs_allres_NoDoublet_slim.RData")

load("/Users/SuganyaSubramanian/MouseCellCycleGenes.rda")

seuObject_slim_nodoub_slim <- CellCycleScoring(seuObject_slim_nodoub_slim,
                                               g2m.features = g2m_genes,
                                               s.features = s_genes)

save(seuObject_slim_nodoub_slim, file ="output_sct_Integrated_Positive/10_SeuratObj_SCT_30pcs_06res_labelled_WithCellCycle.RData")

pdf("output_sct_Integrated_Positive/06_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuObject_slim_nodoub_slim,
        reduction = "umap",
        group.by= "Phase")
dev.off()

posSeurat <- seuObject_slim_nodoub_slim
save(posSeurat, file ="output_sct_Integrated_Positive/10_SeuratObj_SCT_30pcs_06res_labelled_WithCellCycle.RData")
rm(seuObject_slim_nodoub_slim)
