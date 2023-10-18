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
  library(muscat)
  library(Seurat)
  library(clustree)
  #library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(garnett)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(monocle)
})


setwd("/Users/SuganyaSubramanian/Guttridge_Human/")

load("output_sct_Integrated_Negative/08_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered_withHarmony_withCellScoring.RData")

############ FIND MARKERS #############

all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub_withHarmony, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_sct_Integrated_Negative/01_Negative_Presto_Filteredmarkers_padjLT05.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

## top 10 markers for each group
## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_sct_Integrated_Negative/02_Top10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")