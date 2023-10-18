library('rsconnect')
library(ShinyCellPLUS)
library(Seurat)
library(tidyverse)

setwd("/Users/SuganyaSubramanian/Guttridge_Human/")

load("output_Relabel_Positive/12_Positive_FinalShiny.RData")

#https://bioinformatics-musc.shinyapps.io/Pryce-et-al_SingleCell_Human_Cachexia_Positive/

seu_filt <- Pos_Seurat_final

head(seu_filt@meta.data)

seu_filt@meta.data <- seu_filt@meta.data %>% 
  dplyr::select(Cell,Sample,Genotype,seurat_clusters,Phase,nUMI,nGene,pMito,pRibo,
                nCount_RNA,nFeature_RNA,nCount_SCT,nFeature_SCT)

head(seu_filt@meta.data)
DefaultAssay(seu_filt) <- "RNA"
seu_filt <- NormalizeData(object = seu_filt, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)

seu_filt <- FindVariableFeatures(object = seu_filt,
                                 assay = "RNA",
                                 selection.method = "vst")

scConf = createConfig(seu_filt)

makeShinyApp(seu_filt, scConf, gene.mapping = TRUE,
             shiny.title = "Pryce etal SingleCell Human Cachexia Positive",
             markers.all=TRUE, markers.top20=TRUE, de.genes=TRUE, 
             gene.ranks=TRUE, volc.plot=TRUE)

rsconnect::setAccountInfo(name='bioinformatics-musc',
                          token='838A2925138D0715F2D093E909823204',
                          secret='suQLpnyt0hEXRa+LXT9IQxaYNCSIITnH2vBymGPJ')

options(repos = BiocManager::repositories())

rsconnect::deployApp("shinyApp/",
                     appName = 'Pryce-et-al_SingleCell_Human_Cachexia_Positive', 
                     account = 'bioinformatics-musc', 
                     server = 'shinyapps.io')