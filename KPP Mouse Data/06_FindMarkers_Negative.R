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

setwd("/Users/SuganyaSubramanian/BenKpp/")

load("output_sct_Integrated/07_SeuObject_SCT_30pcs_04res_NoDoublet_slim.RData")

############ FIND MARKERS #############

all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub_slim, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_sct_Integrated/Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")







