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
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(Libra)
})

setwd("/Users/SuganyaSubramanian/Guttridge_Human/")

load("output_Relabel/08_SeuratObj_SCT_30pcs_03res_NoDoublet_Reclustered_slim_withCellscoring_labelled_NoEry.RData")

dir.create("output_DGE")

#1 Ctrl Vs Cachetic PDAC
tmp <- subset(seuObject_slim_nodoub_slim_noEry, subset = (ReSample == "Pos_C_PDAC" | ReSample =="Pos_Control"))

table(tmp$ReSample)
#Pos_C_PDAC Pos_Control 
#1073         882

tmp$Sample <- ""
tmp$Sample[tmp$ReSample == "Pos_C_PDAC"] <- "Pos_Exp"
tmp$Sample[tmp$ReSample == "Pos_Control"] <- "Pos_Ctrl"
table(tmp$Sample)
#Pos_Ctrl  Pos_Exp 
#882     1073

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Sample
tmp@meta.data$label <- tmp@meta.data$Sample

DE_Pos_CtrlvsCache <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)
#Neutrophils, Skeletal_Cells,MuSCs_and_Progenitors
#not enough cells, skipping
save(DE_Pos_CtrlvsCache, file = "output_DGE/04_DGE_PosCtrl_vs_PosCachetic.RData")

dge <- DE_Pos_CtrlvsCache %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.0) %>%
  mutate(Direction = case_when(avg_logFC > 0.0 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.0 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/04_DGE_PosCtrl_vs_PosCachetic_Filtered_05_00_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)

#2 Ctrl Vs WS PDAC
tmp <- subset(seuObject_slim_nodoub_slim_noEry, subset = (ReSample == "Pos_WS_PDAC" | ReSample =="Pos_Control"))

table(tmp$ReSample)
#Pos_Control Pos_WS_PDAC 
#882        3725

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$ReSample
tmp@meta.data$label <- tmp@meta.data$ReSample

DE_Pos_CtrlvsWS <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)
# Skeletal_Cells, Pericytes, MuSCs_and_Progenitors
#not enough cells, skipping
save(DE_Pos_CtrlvsWS, file = "output_DGE/05_DGE_PosCtrl_vs_PosWS.RData")

dge <- DE_Pos_CtrlvsWS %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.0) %>%
  mutate(Direction = case_when(avg_logFC > 0.0 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.0 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/05_DGE_PosCtrl_vs_PosWS_Filtered_05_00_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)

#3 Ctrl Vs WS PDAC
tmp <- subset(seuObject_slim_nodoub_slim_noEry, subset = (ReSample == "Pos_WS_PDAC" | ReSample =="Pos_C_PDAC"))

table(tmp$ReSample)
#Pos_C_PDAC Pos_WS_PDAC 
#1073        3725

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$ReSample
tmp@meta.data$label <- tmp@meta.data$ReSample

DE_Pos_CachevsWS <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)
#Neutrophils, Pericytes, MuSCs_and_Progenitors

save(DE_Pos_CachevsWS, file = "output_DGE/06_DGE_PosCache_vs_PosWS.RData")

dge <- DE_Pos_CachevsWS %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.0) %>%
  mutate(Direction = case_when(avg_logFC > 0.0 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.0 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/06_DGE_PosCache_vs_PosWS_Filtered_05_00_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)