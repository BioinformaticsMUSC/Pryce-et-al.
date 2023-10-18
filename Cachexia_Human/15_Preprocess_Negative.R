##prepare seurat for uploading in new shinyplus
setwd("/Users/SuganyaSubramanian/Guttridge_Human/")

load("output_Relabel_Negative/08_Seurat_Negative_Final_Relabeled.RData")

#-----------------
#for adding the gene ranks table output from the AUCell library:
#-----------------

library(Matrix)
library(AUCell)
library(stringr)

Neg_Seurat@misc[['gene_ranks']][['aucell']]$all <- AUCell_buildRankings(Neg_Seurat@assays[['RNA']]@counts)

#--------------------
#for adding the markers through presto after :
#---------------------
##############################
### RUN PRESTO ON CLUSTERS ###
##############################

presto_markers <- presto::wilcoxauc(Neg_Seurat, 'Cell', assay = 'data')
top <- top10 <- presto::top_markers(presto_markers,n = 20,auc_min = 0.5, pval_max = 0.05)
Neg_Seurat@misc[['markers']][['presto']]$all <- presto_markers
Neg_Seurat@misc[['markers']][['presto']]$top_20 <- top

#-------------------------
#for adding diff exp .rds files to the seurat object
#------------------------
library(stringr)

add_diff_exp_file_to_seurat <- function(seurat, file_name, diff_exp_name=NULL) {
  if(!file.exists(file_name)) {
    e<-simpleError(paste0("problem finding file \"", file_name, "\", check the path."))
    stop(e)
  }
  
  split <- tail(strsplit(file_name, "[/]")[[1]], n=1)
  split <- strsplit(split, "[.]")
  
  if( !(tail(split[[1]], n=1) %in% c("rds", "Rdata")) || length(split[[1]]) <= 1) {
    e<-simpleError("file string of type \"file_name.rds\" expected.")
    stop(e)
  }
  if(tail(split[[1]], n=1) == "Rdata") {
    e<-simpleError("file string of type \"file_name.rds\" expected.")
    stop(e)
  }
  # !! call add_diff_exp_obj_to_seurat() here !!
  if(is.null(diff_exp_name)) {
    diff_exp_name = paste0(split[[1]][1:length(split[[1]]) - 1])
  }
  
  data <- readRDS(file_name)
  
  if( !all(class(data) == c('tbl_df', 'tbl', 'data.frame')) ) {
    e<-simpleError("\"tbl_df, tbl, data.frame\" type expected (standard outcome class from Libra's run_de() function: https://github.com/neurorestore/Libra/blob/main/R/run_de.R).")
    stop(e)
  }
  data['de_name'] <- rep(diff_exp_name, times=nrow(data)) # convert as factor?
  if(is.null(seurat@misc$DE_genes$libra$overall)) {
    seurat@misc$DE_genes$libra$overall <- data
  }
  else {
    if( !all(class(seurat@misc$DE_genes$libra$overall) == class(data)) ) {
      e<-simpleError("seurat@misc$DE_genes$libra$overall class type is not the same as data, \"tbl_df, tbl, data.frame\" expected (standard outcome class from Libra's run_de() function: https://github.com/neurorestore/Libra/blob/main/R/run_de.R).")
      stop(e)
    }
    else{
      seurat@misc$DE_genes$libra$overall <- rbind(seurat@misc$DE_genes$libra$overall, data)
    }
  }
  
  return(seurat)
}

load("output_DGE/01_DGE_NegCtrl_vs_NegCachetic.RData")
saveRDS(DE_Neg_CtrlvsCache, file = "output_DGE/DE_Neg_CtrlvsCache.rds")
file1 <- "output_DGE/DE_Neg_CtrlvsCache.rds"

load("output_DGE/02_DGE_NegCtrl_vs_NegWS.RData")
saveRDS(DE_Neg_CtrlvsWS, file = "output_DGE/DE_Neg_CtrlvsWS.rds")
file2 <- "output_DGE/DE_Neg_CtrlvsWS.rds"

load("output_DGE/03_DGE_NegCache_vs_NegWS.RData")
saveRDS(DE_Neg_CachevsWS, file = "output_DGE/DE_Neg_CachevsWS.rds")
file3 <- "output_DGE/DE_Neg_CachevsWS.rds"

Neg_Seurat_final <- add_diff_exp_file_to_seurat(Neg_Seurat,file1,"DE_Neg_CtrlvsCache")
Neg_Seurat_final <- add_diff_exp_file_to_seurat(Neg_Seurat,file2,"DE_Neg_CtrlvsWS")
Neg_Seurat_final <- add_diff_exp_file_to_seurat(Neg_Seurat,file3,"DE_Neg_CachevsWS")

save(Neg_Seurat_final,file = "output_Relabel_Negative/12_Negative_FinalShiny.RData")