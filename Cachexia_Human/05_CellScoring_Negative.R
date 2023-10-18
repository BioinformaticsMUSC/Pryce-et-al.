library(RCurl)
library(AnnotationHub)
library(ensembldb)

load("output_sct_Integrated_Negative/07_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered_withHarmony.RData")

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

seuObject_slim_nodoub_withHarmony <- CellCycleScoring(seuObject_slim_nodoub_withHarmony,
                                                      g2m.features = g2m_genes,
                                                      s.features = s_genes)

pdf("output_sct_Integrated_Negative/13_UMAP_withCellscore.pdf", width = 6, height = 5)
DimPlot(seuObject_slim_nodoub_withHarmony,
        reduction = "umap",
        group.by= "Phase")
dev.off()

save(seuObject_slim_nodoub_withHarmony, file = "output_sct_Integrated_Negative/08_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered_withHarmony_withCellScoring.RData")






