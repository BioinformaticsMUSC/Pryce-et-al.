# Pryce-et-al.
scRNAseq Data Analysis Code
Muscle Inflammation is Regulated by NF-B From Multiple Cells to
Control Distinct States of Wasting in Cancer Cachexia
==========================
This repository contains analysis code for the single cell RNA-seq project carried out by researchers at the [Guttridge lab](https://medicine.musc.edu/departments/pediatrics/clinical-divisions/dcri/faculty/denis-guttridge) and [Berto Lab, MUSC](https://bertolab.org/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: 

## Access the data with an app:

Here a webapp to analyze the data: 
[Pryce etal SingleCell Mouse KPP Negative]([https://bioinformatics-musc.shinyapps.io/Pryce-et-al_SingleCell_Mouse_KPP_Negative/])

## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`input`](input/) | Input/Output data of the initial processing and quality check. | 01_Downstream_Analysis_SCT.r|
| [`output`](output/) | Output data of the initial clustering and integration. | 01_Downstream_Analysis_SCT.r \ 02_Doubletting.r \ 03_Markers_Detection.r|
| [`output_reclust`](output_reclust/) | Output data of the reclustering and subset analyses. | 04_Relabel_And_InitialViz.r \ 05_DGE_Libra_LinearMixed.r \ 07_DGE_visualizations.R \ 09_Database.R \ 10_Motif_Enrichment_scRNA.R|
| [`output_Figure1`](output_Figure1/) | Output for figures and additional visualizations. | 07_DGE_visualizations.R \ 08_GeneOntology.R|
| [`output_SCENIC`](output_SCENIC/) | Output data of the Regulome analyses. | 11_Prepare_Scenic.R \ 12_Analyze_Scenic.R|
| [`output_hdWGCNA`](output_hdWGCNA/) | Output data of the Coexpression analyses. | 13_hdWGCNA.R|
| [`shinyApp`](shinyApp/) | Output of the ShinyApp. | 14_ShinyApp.r|
