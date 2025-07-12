OptiRes: Optimal Resolution Selection for scRNA-seq Analysis
Description
The OptiRes R package is designed to identify the optimal clustering resolution for single-cell RNA sequencing (scRNA-seq) analysis performed with the Seurat package. By leveraging the Silhouette Score algorithm, OptiRes evaluates a range of resolution values (default: 0.01 to 2.00) to determine the resolution that maximizes cluster quality, eliminating the need for arbitrary or experience-based resolution selection. This approach provides a mathematically robust method to enhance the reliability of scRNA-seq clustering results. Additionally, OptiRes includes functionality to visualize silhouette scores and generate colorful dendrograms for cluster relationships, aiding in cell type annotation.
Installation
To install the OptiRes package from GitHub, use the following commands in R:
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install OptiRes from GitHub
devtools::install_github("LingzhangMeng/OptiRes")

# Load OptiRes
library(OptiRes)

Dependencies
The OptiRes package requires the following R packages, which should be installed prior to using OptiRes:
# Install required packages if not already installed
install.packages(c("Seurat", "cluster", "ggplot2", "ggdendro", "dendextend", "circlize"))


Seurat: For scRNA-seq data processing and clustering.
cluster: For computing silhouette scores.
ggplot2: For plotting silhouette scores.
ggdendro: For dendrogram visualization support.
dendextend: For colorful dendrogram plotting.
circlize: For enhanced visualization of cluster relationships.

Tutorial
This tutorial outlines the steps to use OptiRes for determining the optimal clustering resolution and visualizing the results.
Step 1: Preprocess scRNA-seq Data with Seurat
Process your scRNA-seq data using the standard Seurat workflow, which includes the following steps:
# Example Seurat workflow (adjust as needed)
library(Seurat)
seu_obj <- ScaleData(seu_obj)
seu_obj <- RunPCA(seu_obj)
seu_obj <- FindNeighbors(seu_obj)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)  # Temporary resolution
seu_obj <- RunUMAP(seu_obj)

Ensure that UMAP embeddings are generated, as they are required for silhouette score calculations. The initial resolution value in FindClusters can be arbitrary, as OptiRes will evaluate a range of resolutions.
Step 2: Calculate Silhouette Scores Across Resolutions
Use the scSilhouette function to compute silhouette scores for a range of resolution values:
results_df <- scSilhouette(seu_obj, resolutions = seq(0.01, 2.00, 0.01))

Example Output:
Calculating silhouette scores for 200 resolutions...
Computing resolution at 0.01
Silhouette score = 0.3528
Computing resolution at 0.02
Silhouette score = 0.3594
Computing resolution at 0.03
Silhouette score = 0.4341
Computing resolution at 0.04
Silhouette score = 0.4348
...
Computing resolution at 2.00
Silhouette score = 0.3125
Highest Silhouette Score: 0.4348
Optimal Resolution: 0.04

View the results:
View(results_df)

The output is a data frame with two columns: Resolution and SilhouetteScore.
Step 3: Visualize Silhouette Scores
Plot the silhouette scores to identify the optimal resolution:
library(ggplot2)
p <- Plot_res(results_df)
print(p)

This generates a plot with silhouette scores across resolutions, highlighting the optimal resolution in red (e.g., 0.42 in the example dataset).
Step 4: Generate a Colorful Cluster Tree
Create a colorful dendrogram to visualize cluster relationships at the optimal resolution, which aids in cell type annotation:
dendro <- Plot_ColorfulClusterTree(seu_obj, results_df, dims = 1:5)
plot(dendro)

The dims parameter specifies the principal components used for building the cluster tree. Adjust as needed (e.g., dims = 1:3, dims = 1:10) based on your data.
Step 5: Visualize Clusters with UMAP
Visualize the clusters at the optimal resolution using a UMAP plot:
library(Seurat)
DimPlot(seu_obj, reduction = "umap")

This plot displays the cell clusters based on the optimal resolution identified by OptiRes.
Notes

The optimal resolution is automatically selected based on the highest silhouette score, ensuring robust clustering results.
The dims parameter in Plot_ColorfulClusterTree can be adjusted to include more or fewer principal components, depending on the dataset's complexity.
Ensure all dependencies are installed and up-to-date to avoid compatibility issues.

License
This package is licensed under the MIT License. See the LICENSE file for details.
Contact
For questions or issues, please open an issue on the GitHub repository.
