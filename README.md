# Description
The R package OptiRes was developed to figure out the optimal resolution factor frequently used by the Seurat in the process of scRNA-seq analysis,
based on comparison of SilhouetteScore. Previously scientist select resolution from 0.01-2.00, rely on personal experiences, or simply take an arbitray
factor at random. Such a situation usually causes controversial debate over scRNA-seq analysis. While, based on the algorithm Silhoueter Score, this R package - OptiRes,
could provide help choose the optimal resolution factor based on mathematical calculation. Meanwhile, this package will help draw a clustertree for cell type annotations.

# Installation
library(devtools) <br/> 
devtools::install_github("LingzhangMeng/OptiRes")<br/> 
library(OptiRes)<br/> 

# Dependencies
library(Seurat) <br/> 
library(cluster) <br/> 
library(ggplot2) <br/> 
library(ggdendro) <br/> 
library(dendextend) <br/> 
library(circlize) <br/> 

# Tutorial

Step 1: Precess the scRNA-seq data by the standard workflow of R package Seurat, including performing ScaleData(),  RunPCA(), FindNeighbors(), FindClusters() and RunUMAP(). <br/> 
During FindCluster(),take a random resolution factor from 0.01-2.00. Anyway, once the R package OptiRes starts processing analysis, it will calculate Silhouette Score at all resolutons from 0.01-2.00,
and automatically choose the optimal resolution factor with highest Silhouette Score.<br/> <br/> 

Step 2: Calculate Silhoutte Score at different resolutions <br/> 

results_df <- scSilhouette(Seu_obj, resolutions = seq(0.01, 2.00, 0.01)) <br/> 
_________________________________<br/> 
Calculating silhouette scores for 200 resolutions...<br/> 
Calculating silhouette score at resolution: 0.01 <br/> 
  Silhouette score = 0.3528255 <br/> 
Calculating silhouette score at resolution: 0.02 <br/> 
  Silhouette score = 0.3593777 <br/> 
Calculating silhouette score at resolution: 0.03 <br/> 
  Silhouette score = 0.4341156 <br/> 
Calculating silhouette score at resolution: 0.04 <br/> 
  Silhouette score = 0.4347697 <br/> 
Calculating silhouette score at resolution: 0.05 <br/> 
     ..... .....<br/> 
     ..... .....<br/> 
     ..... .....<br/> 
     Calculating silhouette score at resolution: 2.00 <br/> 
_________________________________<br/> <br/> 

View(results_df)
_________________________________<br/> 



_________________________________<br/> <br/> 

Step 3: Plot Silhouette Scores and different resolutons <br/> 

Plot_res(results_df) <br/> 
![Screenshot from 2025-01-02 11-05-36](https://github.com/user-attachments/assets/ed5e1026-e50a-4d11-a445-e4f82a28a4dc)<br/> 

The optimal resoluton factor, 0.42 (for seurat object used in this tutorial) was labelled as red. <br/> <br/> 

Step 4: Draw a colorful cluster tree for checking the relationship between cell types based on genetic profiles,<br/> 
which will help cell type annotations.<br/> 
Plot_ColorfulClusterTree(Seu_obj, results_df) <br/> <br/> 

![Screenshot from 2025-01-02 11-09-31](https://github.com/user-attachments/assets/32f37222-3449-48b6-93d4-88ae57b09ec8)<br/> <br/> 


Step 5: Draw umap at optimal resoluton <br/> 
DimPlot(Seu_obj) <br/> <br/> 
![Screenshot from 2025-01-02 11-15-09](https://github.com/user-attachments/assets/5b5446cf-f82e-4a1e-aa09-a740f5f83d6d)





















































