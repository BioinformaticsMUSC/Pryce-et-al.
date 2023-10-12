##prepare seurat for uploading in new shinyplus
setwd("/Users/SuganyaSubramanian/BenKpp/")

load("Positive/output_Relabel/11_FinalRelabel.RData")

#-----------------
#for adding the gene ranks table output from the AUCell library:
#-----------------

library(Matrix)
library(AUCell)
library(stringr)

negSeurat@misc[['gene_ranks']][['aucell']]$all <- AUCell_buildRankings(negSeurat@assays[['RNA']]@counts)

#--------------------
#for adding the markers through presto after :
#---------------------
##############################
### RUN PRESTO ON CLUSTERS ###
##############################

presto_markers <- presto::wilcoxauc(negSeurat, 'Cell', assay = 'data')
top <- top10 <- presto::top_markers(presto_markers,n = 20,auc_min = 0.5, pval_max = 0.05)
negSeurat@misc[['markers']][['presto']]$all <- presto_markers
negSeurat@misc[['markers']][['presto']]$top_20 <- top

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

load("output_DGE/01_DGE_NegCtrl_vs_NegKpp.RData")
saveRDS(DE_Neg_CtrlvsKpp, file = "output_DGE/DE_Neg_CtrlvsKpp.rds")
file1 <- "output_DGE/DE_Neg_CtrlvsKpp.rds"

negSeurat_final <- add_diff_exp_file_to_seurat(negSeurat,file1,"DE_Neg_CtrlvsKpp")

save(negSeurat_final,file = "Positive/output_Relabel/12_Positive_FinalShiny.RData")
