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

load("output_Relabel_Negative/08_Seurat_Negative_Final_Relabeled.RData")

dir.create("output_DGE")

table(Neg_Seurat$Sample)
#Neg_C_PDAC Neg_Control Neg_WS_PDAC 
#       4368        5463        1696 

#1 Ctrl Vs Cachetic PDAC
tmp <- subset(Neg_Seurat, subset = (Sample == "Neg_C_PDAC" | Sample =="Neg_Control"))

table(tmp$Sample)
#Neg_C_PDAC Neg_Control 
#4368        5463

tmp$Sample[tmp$Sample == "Neg_C_PDAC"] <- "Neg_Exp"
tmp$Sample[tmp$Sample == "Neg_Control"] <- "Neg_Ctrl"
table(tmp$Sample)

#Neg_Ctrl  Neg_Exp 
#5529     4404

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Sample
tmp@meta.data$label <- tmp@meta.data$Sample

DE_Neg_CtrlvsCache <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

save(DE_Neg_CtrlvsCache, file = "output_DGE/01_DGE_NegCtrl_vs_NegCachetic.RData")

dge <- DE_Neg_CtrlvsCache %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.2) %>%
  mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/01_DGE_NegCtrl_vs_NegCachetic_Filtered_05_02_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)

#2 Ctrl Vs WS PDAC
tmp <- subset(Neg_Seurat, subset = (Sample == "Neg_WS_PDAC" | Sample =="Neg_Control"))

table(tmp$Sample)
#Neg_Control Neg_WS_PDAC 
#5390        1290

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Sample
tmp@meta.data$label <- tmp@meta.data$Sample

DE_Neg_CtrlvsWS <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

save(DE_Neg_CtrlvsWS, file = "output_DGE/02_DGE_NegCtrl_vs_NegWS.RData")

dge <- DE_Neg_CtrlvsWS %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.2) %>%
  mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/02_DGE_NegCtrl_vs_NegWS_Filtered_05_02_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)

#3 Ctrl Vs WS PDAC
tmp <- subset(Neg_Seurat, subset = (Sample == "Neg_WS_PDAC" | Sample =="Neg_C_PDAC"))

table(tmp$Sample)
#Neg_C_PDAC Neg_WS_PDAC 
#4368        1696

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Sample
tmp@meta.data$label <- tmp@meta.data$Sample

DE_Neg_CachevsWS <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

save(DE_Neg_CachevsWS, file = "output_DGE/03_DGE_NegCache_vs_NegWS.RData")

dge <- DE_Neg_CachevsWS %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.2) %>%
  mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/03_DGE_NegCache_vs_NegWS_Filtered_05_02_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)