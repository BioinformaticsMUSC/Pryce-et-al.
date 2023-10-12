# Pryce-et-al.
scRNAseq Data Analysis Code
Muscle Inflammation is Regulated by NF-kB From Multiple Cells to
Control Distinct States of Wasting in Cancer Cachexia
==========================
This repository contains analysis code for the single cell RNA-seq project carried out by researchers at the [Guttridge lab, MUSC](https://medicine.musc.edu/departments/pediatrics/clinical-divisions/dcri/faculty/denis-guttridge) and [Berto Lab, MUSC](https://bertolab.org/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: 

## Access the data with an app:

Here a webapp to analyze the data:

[Pryce etal SingleCell Mouse KPP Negative](https://bioinformatics-musc.shinyapps.io/Pryce-et-al_SingleCell_Mouse_KPP_Negative/)

[Pryce etal SingleCell Mouse KPP Positive](https://bioinformatics-musc.shinyapps.io/Pryce-et-al_SingleCell_Mouse_KPP_Positive/)

## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`KPP_MouseData_Analysis`](KPP_MouseData_Analysis/) | Code and Figuers for Mouse KPP Data | |
| [`output_Positive`](KPP_MouseData_Analysis/output_sct_Integrated_Positive) | Output data of the initial clustering and integration of Positive Dataset - CtrlPos and KppPos | 01_Seurat_Clustering_Positive.R \ 03_Doubletting_Positive.R \ 05_FindMarkers_Positive.R \ 07_CellCycleScoring.R|
| [`output_relabel_Positive`](KPP_MouseData_Analysis/output_Relabel_Positive/) | Output data of Cell annotation | 08_Relable_Positive.R| 
| [`output_Negative`](KPP_MouseData_Analysis/output_sct_Integrated_Negative) | Output data of the initial clustering and integration of Negative Dataset - CtrlNeg and KppNeg | 02_Seurat_Clustering_Negative.R \ 04_Doubletting_Negative.R \ 06_FindMarkers_Negative.R \ 07_CellCycleScoring.R|
| [`output_relabel_Negative`](KPP_MouseData_Analysis/output_Relabel_Negative/) | Output data of Cell annotation | 09_Relable_Negative.R|
| [`output_DGE`](KPP_MouseData_Analysis/output_DGE/) | Output of Differential gene expression for Positive and Negative data | 10_DGE.R |
| [`output_Figure_Positive`](KPP_MouseData_Analysis/output_Figure_Positive/) | Output data of final figures and FAPs subset analysis. | 11_Final_Figuers_Positive.R|
| [`output_Figure_Negative`](KPP_MouseData_Analysis/output_Figure_Negative/) | Output data of final figures and Macrophage subset analysis. | 11_Final_Figuers_Negative.R|
| [`shinyApp`](KPP_MouseData_Analysis/ShinyApp/) | Output of the ShinyApp. | 13_PreprocessForShiny_Positive.R \ 14_ShinyApp_Positive.R \ 15_PreprocessForShiny_Negative.R \ 16_ShinyApp_Negative.R |
