# OptiRes: Optimal Resolution Selection for scRNA-seq Analysis

## Description

The `OptiRes` R package is designed to identify the optimal clustering resolution for single-cell RNA sequencing (scRNA-seq) analysis performed with the `Seurat` package. By leveraging the Silhouette Score algorithm, `OptiRes` evaluates a range of resolution values (default: 0.01 to 2.00) to determine the resolution that maximizes cluster quality, eliminating the need for arbitrary or experience-based resolution selection. This approach provides a mathematically robust method to enhance the reliability of scRNA-seq clustering results. Additionally, `OptiRes` includes functionality to visualize silhouette scores and generate colorful dendrograms for cluster relationships, aiding in cell type annotation. By the way, this package is also suitable for integrated seurat objects.

## Installation

**To install the `OptiRes` package from GitHub, use the following commands in R.** 

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install OptiRes from GitHub
devtools::install_github("LingzhangMeng/OptiRes")

# Load OptiRes
library(OptiRes)
```

## Dependencies

The `OptiRes` package requires the following R packages, which should be installed prior to using `OptiRes`. Copy the code below to install dependencies:

```R
# Install required packages if not already installed
install.packages(c("Seurat", "cluster", "ggplot2", "ggdendro", "dendextend", "circlize"))

# Load dependencies
library(Seurat)
library(cluster)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(circlize)
library(crayon)
```

- `Seurat`: For scRNA-seq data processing and clustering.
- `cluster`: For computing silhouette scores.
- `ggplot2`: For plotting silhouette scores.
- `ggdendro`: For dendrogram visualization support.
- `dendextend`: For colorful dendrogram plotting.
- `circlize`: For enhanced visualization of cluster relationships.

## Tutorial

This tutorial outlines the steps to use `OptiRes` for determining the optimal clustering resolution and visualizing the results. You can copy each code snippet by clicking the copy icon that appears when hovering over the code block on GitHub.

### Step 1: Preprocess scRNA-seq Data with Seurat

**Process your scRNA-seq data using the standard `Seurat` workflow:**

```R
# Example Seurat workflow (adjust as needed)
library(Seurat)
seu_obj <- ScaleData(seu_obj)
seu_obj <- RunPCA(seu_obj)
seu_obj <- FindNeighbors(seu_obj)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)  # Temporary resolution
seu_obj <- RunUMAP(seu_obj)
```

Ensure that UMAP embeddings are generated, as they are required for silhouette score calculations. The initial resolution value in `FindClusters` can be arbitrary, as `OptiRes` will evaluate a range of resolutions.

### Step 2: Calculate Silhouette Scores Across Resolutions

**Use the `scSilhouette` function to compute silhouette scores for a range of resolution values:**

```R
results_df <- scSilhouette(seu_obj, resolutions = seq(0.01, 2.00, 0.01))
```

**Example Output:**

```
Calculating silhouette scores for 200 resolutions...
Computing resolution at 0.01
Silhouette score = 0.3496
Computing resolution at 0.02
Silhouette score = 0.4179
Computing resolution at 0.03
Silhouette score = 0.4175
Computing resolution at 0.04
Silhouette score = 0.4522
Computing resolution at 0.05
Silhouette score = 0.4787
Computing resolution at 0.06
Silhouette score = 0.4862
Computing resolution at 0.07
Silhouette score = 0.4864
Computing resolution at 0.08
Silhouette score = 0.4819
Computing resolution at 0.09
Silhouette score = 0.4726
Computing resolution at 0.10
Silhouette score = 0.3838
Computing resolution at 0.11
...
Computing resolution at 2.00
Silhouette score = 0.3125
Highest Silhouette Score: 0.4864
Optimal Resolution: 0.07
```

View the results:

```R
View(results_df)
```

The output is a data frame with two columns: `Resolution` and `SilhouetteScore`.

```R
# Abstract the value of the optimal resolution
optimal_resolution = results_df$Resolution[which.max(results_df$SilhouetteScore)]
cat("Optimal Resolution =", red(optimal_resolution))
```
**Example Output:**

```
Optimal Resolution = 0.07
```

### Step 3: Visualize Silhouette Scores

**Plot the silhouette scores to identify the optimal resolution:**

```R
p <- Plot_res(results_df)
print(p)
```
**Example Output:**
<img width="931" height="704" alt="Screenshot 2025-07-12 092636" src="https://github.com/user-attachments/assets/c23bba10-f352-4c8a-a332-5b6b3d5215af" />

This generates a plot with silhouette scores across resolutions, highlighting the optimal resolution in red (e.g., 0.07 in the example dataset).

### Step 4: Generate a Colorful Cluster Tree

**Create a colorful dendrogram to visualize cluster relationships at the optimal resolution, which aids in cell type annotation:**

```R
dendro <- Plot_ColorfulClusterTree(seu_obj, results_df, dims = 1:5)
plot(dendro)
```

**Example Output:**

<img width="393" height="363" alt="Weixin Image_20250713103018" src="https://github.com/user-attachments/assets/bae5b73b-5fd6-4818-a20a-aefdcd7e7ea4" />


The `dims` parameter specifies the principal components used for building the cluster tree. Adjust as needed (e.g., `dims = 1:3`, `dims = 1:10`) based on your data.

### Step 5: Visualize Clusters with UMAP

**Visualize the clusters at the optimal resolution using a UMAP plot:**

```R
# Re-cluster on the optimal resolution factor
Cell.integrated <- FindClusters(seu_obj, pc.use = 1:10, resolution = optimal_resolution, 
                                                 group.singletons = TRUE, verbose = 0, save.SNN = T)
DimPlot(seu_obj, reduction = "umap")
```

This plot displays the cell clusters based on the optimal resolution identified by `OptiRes`.

## Notes

- The optimal resolution is automatically selected based on the highest silhouette score, ensuring robust clustering results.
- The `dims` parameter in `Plot_ColorfulClusterTree` can be adjusted to include more or fewer principal components, depending on the dataset's complexity.
- Ensure all dependencies are installed and up-to-date to avoid compatibility issues.
- **To copy code snippets, hover over the code block on GitHub and click the copy icon that appears in the top-right corner.**

## License

This package is licensed under the MIT License. See the LICENSE file for details.

## Contact

For questions or issues, please open an issue on the [GitHub repository](https://github.com/LingzhangMeng/OptiRes).
