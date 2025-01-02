# Description
The R package OptiRes was developed to figure out the optimal resolution factor frequently used by the Seurat in the process of scRNA-seq analysis,
based on comparison of SilhouetteScore. Previously scientist select resolution from 0.01-2.00, rely on personal experiences, or simply take an arbitray
factor at random. Such a situation usually causes controversial debate over scRNA-seq analysis. While, based on the algorithm Silhoueter Score, this R package - OptiRes,
could provide help choose the optimal resolution factor based on mathematical calculation.

# Installation
library(devtools)
devtools::install_github("LingzhangMeng/OptiRes")
library(OptiRes)
