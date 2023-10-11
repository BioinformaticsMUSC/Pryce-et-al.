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
  library(Libra) #devtools::install_github("neurorestore/Libra")
})

setwd("/Users/SuganyaSubramanian/BenKpp/")
dir.create("output_DGE")

#### Negative data #####
load("Negative/output_Relabel/09_SeuratObj_SCT_30pcs_04res_labelled.RData")
load("Negative/output_Relabel/11_FinalRelabel.RData")

table(negSeurat$Genotype)

tmp <- negSeurat
tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Genotype
tmp@meta.data$label <- tmp@meta.data$Genotype

DE_Neg_CtrlvsKpp <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

save(DE_Neg_CtrlvsKpp, file = "output_DGE/01_DGE_NegCtrl_vs_NegKpp.RData")

dge <- DE_Neg_CtrlvsKpp %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.0) %>%
  mutate(Direction = case_when(avg_logFC > 0.0 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.0 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/01_DGE_NegCtrl_vs_NegKpp_Filtered_05_00_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)

##### Positive #####
#load("Positive/output_Relabel/09_SeuratObj_SCT_30pcs_06res_labelled.RData")
load("Positive/output_Relabel/11_FinalRelabel.RData")

table(posSeurat$Genotype)

tmp <- posSeurat
tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Genotype
tmp@meta.data$label <- tmp@meta.data$Genotype

DE_Pos_CtrlvsKpp <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

save(DE_Pos_CtrlvsKpp, file = "output_DGE/02_DGE_PosCtrl_vs_PosKpp.RData")

dge <- DE_Pos_CtrlvsKpp %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.0) %>%
  mutate(Direction = case_when(avg_logFC > 0.0 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.0 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(dge, 
                     file = "output_DGE/02_DGE_PosCtrl_vs_PosKpp_Filtered_05_00_withDirection.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)