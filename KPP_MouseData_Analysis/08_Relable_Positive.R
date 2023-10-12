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
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(RColorBrewer)
  library(MetBrewer)
  library(scales)
  library(ComplexHeatmap)
  library(randomcoloR)
  library(png)#to include pic in ggplot
  library(patchwork)
})

setwd("/Users/SuganyaSubramanian/BenKpp/")

load("output_Relabel_Positive/10_SeuratObj_SCT_30pcs_06res_labelled_WithCellCycle.RData")

###2. Load Labels data
Labels <- read.table("output_sct_Integrated_Positive/New_Labels.txt",header=T,sep="\t")


#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

current.cluster.ids
new.cluster.ids

#Idents(posSeurat) <- "seurat_clusters"
#table(Idents(posSeurat))

##Step 3: Rename cluster with new Labels
posSeurat@active.ident <- plyr::mapvalues(x = posSeurat@active.ident, 
                                          from = current.cluster.ids, 
                                          to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
posSeurat@meta.data$Cell <- posSeurat@active.ident 

##Save new Seurat data
save(posSeurat,file = "output_Relabel_Positive/11_FinalRelabel.RData")

pdf("output_Relabel_Positive/01_UMAP_Labelled.pdf", width = 6, height = 7)
DimPlot(object = posSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("output_Relabel_Positive/02_UMAP_Labelled_Groupby.pdf", width = 12, height = 7)
p1 <- DimPlot(object = posSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = posSeurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel_Positive/03_UMAP_Labelled_withclusters.pdf", width = 12, height = 7)
Idents(posSeurat) <- "seurat_clusters"
p1 <- DimPlot(object = posSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(posSeurat) <- "Cell"
p2 <- DimPlot(object = posSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel_Positive/04_UMAP_Labelled_Splitby_Genotype.pdf", width = 12, height = 7)
DimPlot(object = posSeurat, reduction = "umap", label = TRUE,ncol = 2, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Faps","Skeletal","Schwann","Pericytes","Endothelial",
           "B_Cells","T_Cells","Macrophage","Neutrophil",
           "Basophil","Myeloid_DC","Erythrocytes","Platelets") 

Idents(posSeurat) <- factor(Idents(posSeurat), levels = order) # reorder the factor based on cell type


markers <- c("Gsn",#Faps "Col1a1",
             "Tnnc2",#skeletal, MuSC-"Pax7",
             "Pou3f1",#Schwann
             "Acta2",#pericytes
             "Pecam1", #endo, Epi-"Lgals7","Lor",
             "Ms4a1",#Bcells
             "Cd3d",#Tcells
             #"Arhgap15",#Tcells(?)
             "C1qa",#Macrophage
             "S100a9",#Neutrophil
             "Cma1",#Basophil
             "Fscn1",#MyeloidDC
             "Hbb-bs","Pf4") #Erythrocytes, Platelets

pdf("output_Relabel_Positive/07_DotPlot_Celltype_Markers.pdf", width = 10, height = 6)
DotPlot(posSeurat, features = markers) + rotate_x_text(angle = 45)
dev.off()


































##### Save Alternative PLOTS ######
### 1 UMAP ####
umap <- as.data.frame(Embeddings(posSeurat, reduction = "umap"))
meta <- as.data.frame(posSeurat@meta.data)

df <- cbind(umap,meta)%>% 
  dplyr::group_by(Cell) %>% 
  dplyr::mutate(N = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Genotype) %>%
  as.data.frame()

label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  dplyr::group_by(Cell) %>% 
  dplyr::summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  dplyr::left_join(label) %>%
  as.data.frame() 

#choose No. depending on the total no. of clusters
col <- randomcoloR::distinctColorPalette(24)

write.table(col, "output_Relabel/Colors_Used.txt", sep="\t",quote=F)

#Manually edited the colors
#col <- data.table::fread("output_sct_0.6/Colors_Used_Initial_Data_edited.txt")
#col <- as.character(col$V2)
#col

pdf("output_Relabel_Positive/05_UMAP_Labelled_Alternative.pdf", width = 8, height = 7)
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = Cell),size=0.2) +
  ggrepel::geom_text_repel(data = label_2, aes(label = label),
                           color = "black",
                           #fontface = 'bold',
                           segment.colour = "grey60",
                           box.padding = unit(0.25, "lines"),
                           point.padding = unit(0.5, "lines"),
                           nudge_x = .15,
                           nudge_y = 1,
                           size = 6) + 
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  scale_colour_manual(values = col)+
  theme_void() +
  theme(legend.position="none")
dev.off()

### 2 Cell proportion ####

color <- c("#17154f","#b0799a","#e69b00","#355828","#4f2400","#b62684")
#"#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

prop_cell <-  posSeurat@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Cell) %>%
  dplyr::group_by(Cell,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Cell) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Cell) %>%
  ggbarplot("Cell", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_Relabel/06_CellProportion_By_Genotype.pdf", plot = prop_cell, width = 16, height = 18, units = "in", dpi = 150)


# 03 Cell Count
cnt <- table(posSeurat$Cell, posSeurat$Genotype)
openxlsx::write.xlsx(cnt,file = "output_Relabel/Cell_Count.xlsx", colNames = TRUE, borders = "columns")

#Dotplot

order <- c("Y143_Cells_1","Y143_Cells_2","Y143_Cells_3",
           "LLC_Cells_1","LLC_Cells_2","Tumor_Cells_1","Tumor_Cells_2",
           "M1_Macrophages_1","M1_Macrophages_2",
           "M2_Macrophages_1","M2_Macrophages_2",
           "B_Cells","T_Cells","Neutrophils","Endothelial_Cells",
           "Keratin_Expressing_Stroma_Cells_1","Keratin_Expressing_Stroma_Cells_2",
           "Cancer_Associated_Fibroblasts","Muscle_Cells","UnDefined1","UnDefined2")


markers <- c("Mgp","Spp1","Clu",
             "Rhox5","Cdh18","Rpl37","Dst",
             "Lyz2","Cd74","C1qa","Clec4a2",
             "Cd79a","Cd3e","Itk","Nkg7","S100a9","Pecam1",
             "Lgals7","Dsc3",
             "Col1a2","Des")

group <- c("Spp1","Rhox5","Rpl37","Dst",
           "Lyz2","Cd74","C1qa","Clec4a2",
           "Cd79a","Itk","S100a9","Pecam1",
           "Dsc3",
           "Col1a2","Des")

Idents(posSeurat) <- factor(Idents(posSeurat), levels = order) # reorder the factor based on cell type

pdf("output_Relabel_Positive/07_DotPlot_Celltype_Markers.pdf", width = 10, height = 6)
DotPlot(posSeurat, features = markers) + rotate_x_text(angle = 45)
dev.off()

pdf("output_Relabel_Positive/07_DotPlot_Celltype_MarkersGroup.pdf", width = 10, height = 6)
DotPlot(posSeurat, features = group) + rotate_x_text(angle = 45)
dev.off()

#Rhox5 Gene expression Analysis
Rhox5 = GetAssayData(object=posSeurat,assay="RNA",slot="data")["Rhox5",]

posSeurat@meta.data$Rhox5 <- as.numeric(log2(Rhox5+1))

umap <- as.data.frame(Embeddings(posSeurat, reduction = "umap"))
meta <- as.data.frame(posSeurat@meta.data)

# filter Rhox5 only
df <- cbind(umap,meta)%>%
  filter(Rhox5 > 0)

grouped_data <- table(df$Cell, df$Genotype) %>%
  as.matrix() %>%
  as.data.frame() %>%
  arrange(Var1,Var2) %>%
  dplyr::rename(Cell = Var1, Genotype = Var2, Cells = Freq)  %>%
  mutate(Cell=fct_rev(fct_relevel(Cell,order))) %>%
  as.data.frame() %>%
  droplevels()

openxlsx::write.xlsx(grouped_data,file = "output_Relabel/Rhox5_data.xlsx", colNames = TRUE, borders = "columns")

pdf("output_Relabel_Positive/08_Barplot_NumberOfCells_Rhox5.pdf", width = 8, height = 10)
ggbarplot(grouped_data, "Cell", "Cells",
          fill = "Genotype", color = "Genotype", palette = color,
          label = FALSE, lab.col = "white", lab.pos = "in",
          rotate = TRUE) +
  #scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="top") +
  xlab("")
dev.off()

pdf("output_Relabel_Positive/08_Barplot_NumberOfCells_Rhox5_Hori.pdf", width = 10, height = 6)
ggbarplot(grouped_data, "Cell", "Cells",
          fill = "Genotype", color = "Genotype", palette = color,
          label = FALSE, lab.col = "white", lab.pos = "in"#,rotate = TRUE
) +
  #scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")
dev.off()

#Feature plot
df <- cbind(umap,meta)

pdf("output_Relabel_Positive/09_FeaturePlot_NumberOfCells_Rhox5.pdf", width = 8, height = 7)
ggplot(df, aes(x=UMAP_1, y=UMAP_2,color = Rhox5)) +
  #geom_point() +
  ggrastr::geom_point_rast(aes(colour = Rhox5),size=0.2) +
  scale_colour_gradient2(low = "gray60", high = "darkblue", mid="gray60",na.value = NA) +##FDE725FF #"#440154FF"
  theme_void() +
  theme(legend.position="none")
dev.off()

#Presto Markers by Cell type
all_markers_clustID <- presto::wilcoxauc(posSeurat, 'Cell', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

#openxlsx::write.xlsx(all_markers_clustID, file = "output_Relabel/Presto_Markers_All_CellType.xlsx", colNames = TRUE, borders = "columns")

all_markers_clustID.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05)

openxlsx::write.xlsx(all_markers_clustID.Sign, file = "output_Relabel/Presto_Markers_filtered_padjlt05_CellType.xlsx", colNames = TRUE, borders = "columns")

#Can try
#https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)
#DefaultAssay(posSeurat)  <- "RNA" 
#VariableFeaturePlot_scCustom(seurat_object = posSeurat, 
#                             num_features = 20, repel = TRUE,assay = "RNA",
#                             y_axis_log = TRUE)
#Error: Unable to find highly variable feature information for method 'sct'

# Set color palette
pal <- viridis(n = 10, option = "D")
pdf("output_Relabel_Positive/09_FeaturePlot_NumberOfCells_Rhox5_scCustom.pdf", width = 8, height = 7)
FeaturePlot_scCustom(seurat_object = posSeurat, 
                     features = "Rhox5", colors_use = pal,raster = "TRUE")
dev.off()

pdf("output_Relabel_Positive/10_DensityPlot_Rhox5.pdf", width = 8, height = 7)
Plot_Density_Custom(seurat_object = posSeurat, 
                    features = "Rhox5")
#custom_palette = PurpleAndYellow())
dev.off()
##### Not making these
### Heatmaps by cell type
posSeurat.markers <- FindAllMarkers(object = posSeurat, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25, random.seed = 1, return.thresh = 0.01)

write.table(data.frame(posSeurat.markers), file = "output_Relabel/integrated_markers_per_cluster.txt", 
            col.names = T, row.names = F, sep = "\t", quote = F)

top.markers <- as.data.frame(posSeurat.markers %>% group_by(cluster) %>% top_n(30, avg_log2FC))

write.table(top.markers, file = "output_Relabel/integrated_markers_per_cluster.top30.txt", 
            col.names = T, row.names = F, sep = "\t", quote = F)

posSeurat <- FindVariableFeatures(object = posSeurat,
                                  assay = "RNA",
                                  selection.method = "vst")

posSeurat <- ScaleData(posSeurat,vars.to.regress = c("nUMI","pMito"))

posSeurat <- RunPCA(object = posSeurat, 
                    features=NULL, 
                    weight.by.var = TRUE, 
                    ndims.print = 1:5, 
                    nfeatures.print = 30, 
                    npcs = 30, 
                    reduction.name = "pca")

#P2 <- DoHeatmap(posSeurat, features = unique(top.markers_5$gene), 
#                cells = sample(colnames(posSeurat),5000), group.by = "ident", 
#                group.bar = TRUE, disp.min = -2.5,disp.max = NULL,
#                slot = "scale.data", 
#                assay = NULL, label = TRUE, size = 2, hjust = 0, angle = 45, 
#                raster = FALSE, draw.lines = TRUE, lines.width = NULL, 
#                #group.colors = c("salmon4","darkslategray4", "salmon2", "palevioletred1", "darkslateblue", "tomato", 
#                #                 "darkorange", "khaki", "red", "goldenrod1", "red4", "honeydew3"), 
#                group.bar.height = 0.02,
#                combine = TRUE)
#P2 <- P2 + theme(aspect.ratio = 4/1)
#P2 <- P2 + theme(text = element_text(size = 6))
#pdf("output_Relabel_Positive/08_HeatMap_Celltype_Markers.pdf", width = 8, height = 10)
#P2
#dev.off()

top.markers_5 <- as.data.frame(posSeurat.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC))

P2 <- DoHeatmap(posSeurat, features = unique(top.markers_5$gene), 
                cells = sample(colnames(posSeurat),5000), group.by = "ident", 
                group.bar = TRUE, disp.min = -2.5,disp.max = NULL,
                slot = "scale.data", 
                assay = NULL, label = TRUE, size = 2, hjust = 0, angle = 45, 
                raster = FALSE, draw.lines = TRUE, lines.width = NULL, 
                #group.colors = c("salmon4","darkslategray4", "salmon2", "palevioletred1", "darkslateblue", "tomato", 
                #                 "darkorange", "khaki", "red", "goldenrod1", "red4", "honeydew3"), 
                group.bar.height = 0.02,
                combine = TRUE)

P2 <- P2 + theme(aspect.ratio = 4/1)
P2 <- P2 + theme(text = element_text(size = 6))
pdf("output_Relabel_Positive/08_HeatMap_Celltype_Markers_top5.pdf", width = 8, height = 15)
P2
dev.off()

