suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(sctransform)
  library(hdf5r)
  library(ggrastr)
  library(clustree)
  library(cowplot)
})

setwd("/Users/SuganyaSubramanian/BenKpp/")


CtrlNeg <- Read10X_h5("cellbender_out/Ctrl_Neg_filtered.h5", use.names = TRUE, unique.features = TRUE)
KppNeg <- Read10X_h5("cellbender_out/Kpp_Neg_filtered.h5", use.names = TRUE, unique.features = TRUE)

CtrlNeg_obj <- CreateSeuratObject(counts = CtrlNeg,min.features = 100)
CtrlNeg_obj

KppNeg_obj <- CreateSeuratObject(counts = KppNeg,min.features = 100)
KppNeg_obj


seuObject <- merge(CtrlNeg_obj,
                   y = c(KppNeg_obj),
                   add.cell.ids = c("CtrlNeg","KppNeg"),
                   project = "A9F1")

seuObject[["pMito"]] <- PercentageFeatureSet(seuObject, pattern = "^mt-")
seuObject[["pRibo"]] <- PercentageFeatureSet(seuObject,pattern = "^Mrp[sl]")

seuObject@meta.data <- seuObject@meta.data %>%
  rownames_to_column("Cell") %>%
  mutate(Genotype = sapply(X = strsplit(colnames(seuObject), split = "_"), FUN = "[", 1)) %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
  column_to_rownames("Cell")

head(seuObject@meta.data)

##########STEP6############# Save Seurat data

dir.create("output_sct_Integrated_Negative")

save(seuObject,file="output_sct_Integrated_Negative/01_SeuratObj_Unfilt.RData")

##########STEP7#############check data with plots
#QC plot 1
pdf("output_sct_Integrated_Negative/01_Quality_Control_plots.pdf", width=6,height=4)
feats <- c("nUMI", "nGene", "pMito")
VlnPlot(seuObject, group.by = "Genotype", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
dev.off()

pdf("output_sct_Integrated_Negative/01_Quality_Control_plots_log.pdf", width=6,height=4)
feats <- c("nUMI", "nGene", "pMito")
VlnPlot(seuObject, group.by = "Genotype", features = feats,log = TRUE, pt.size = 0, ncol = 3) + 
  NoLegend()
dev.off()

##-------------------------------------------------------
## DATA FILTERING - Step 9
##-------------------------------------------------------

#####################
## Remove MT genes ##
#####################
mito_filtered <- seuObject@assays$RNA@counts[-grep("^mt-",rownames(seuObject@assays$RNA@counts)),]

# Initialize the Seurat object with the raw (non-normalized data).
seuObject_final <- CreateSeuratObject(counts = mito_filtered, project = "A9F1")

## Add pMito info from meta data for all cells before filtering
metaAll <- as.data.frame(seuObject@meta.data)
seuObject_final <- AddMetaData(object = seuObject_final, metadata = as.data.frame(seuObject@meta.data))
seuObject_final@meta.data$nCount_RNA <- NULL
seuObject_final@meta.data$nFeature_RNA <- NULL

save(seuObject_final,file="output_sct_Integrated_Negative/02_SeuratObj_final.RData")

seuObject_filt <- subset(x = seuObject_final, subset = nUMI < 10000 & pMito < 20)

pdf("output_sct_Integrated_Negative/01_Quality_Control_plots_filtered.pdf", width=6,height=4)
feats <- c("nUMI", "nGene", "pMito")
VlnPlot(seuObject_filt, group.by = "Genotype", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
dev.off()

save(seuObject_filt,file="output_sct_Integrated_Negative/02_SeuratObj_Filt.RData")

#Gene Expression Barplot
df <- table(seuObject_filt$Genotype) %>% as.data.frame()

pdf("output_sct_Integrated_Negative/06_GeneExp_barplot.pdf",width = 5,height = 5)
ggplot(df,aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  #scale_fill_manual(values=c("orange2", "green3","red1")) +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  ggtitle("Gene Expression by Sample") +
  xlab("Sample") +
  ylab("Gene Expression") +
  theme(legend.position="none")
dev.off()

##-------------------------------------------------------
## DATA INTEGRATION
##-------------------------------------------------------

seuObject_split <- SplitObject(seuObject_filt, split.by = "Genotype")

seuObject_split <- seuObject_split[c("CtrlNeg","KppNeg")]

for (i in 1:length(seuObject_split)) {
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]], 
                                      vars.to.regress = c("nUMI","pMito","pRibo"),
                                      verbose = FALSE)
}

integ_features <- SelectIntegrationFeatures(object.list = seuObject_split, 
                                            nfeatures = 4000) 


seuObject_split <- PrepSCTIntegration(object.list = seuObject_split, 
                                      anchor.features = integ_features)


integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

save(integ_anchors, file = "output_sct_Integrated_Negative/tmp3.RData")

seuObject_integrated <- IntegrateData(
  anchorset = integ_anchors,
  new.assay.name = "integrated",
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 100,
  sd.weight = 1,
  eps = 0.5,
  verbose = TRUE
)
save(seuObject_integrated, file = "output_sct_Integrated_Negative/tmp4.RData")

DefaultAssay(seuObject_integrated) <- "integrated"

seuObject_integrated <- RunPCA(object = seuObject_integrated, 
                               features=NULL, 
                               weight.by.var = TRUE, 
                               ndims.print = 1:5, 
                               nfeatures.print = 30, 
                               npcs = 30, 
                               reduction.name = "pca")

seuObject_integrated <- FindNeighbors(object = seuObject_integrated, 
                                      reduction = "pca", 
                                      dims = 1:30, 
                                      nn.eps = 0.5)

seuObject_integrated <- FindClusters(object = seuObject_integrated, 
                                     resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2), 
                                     algorithm = 1,
                                     n.iter = 1000)

save(seuObject_integrated, file = "output_sct_Integrated_Negative/03_SeuObj_Integrated_allres.RData")

pdf("output_sct_Integrated_Negative/07_Data_Integrated_Clustree.pdf", width = 12, height = 6)
clustree(seuObject_integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Select resolution and run UMAP
Idents(object = seuObject_integrated) <- "integrated_snn_res.0.5"

seuObject_integrated <- RunUMAP(object = seuObject_integrated, 
                                reduction = "pca", 
                                dims = 1:30)

# Select the RNA counts slot to be the default assay
DefaultAssay(seuObject_integrated) <- "RNA"

seuObject_integrated <- NormalizeData(object = seuObject_integrated, 
                                      normalization.method = "LogNormalize", 
                                      scale.factor = 10000)

save(seuObject_integrated, file = "output_sct_Integrated_Negative/04_SeuratObj_SCT_30pcs_05res.RData")


pdf("output_sct_Integrated_Negative/08_UMAP_Data_Integrated.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

seuObject_integrated_slim <- DietSeurat(seuObject_integrated, 
                                        counts = TRUE, 
                                        data = TRUE, 
                                        scale.data = FALSE,
                                        assays="RNA",
                                        dimreducs = c("pca","umap"))

save(seuObject_integrated_slim, file = "output_sct_Integrated_Negative/05_SeuratObj_SCT_30pcs_05res_Slimmed.RData")
