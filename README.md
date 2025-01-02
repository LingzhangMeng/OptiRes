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
