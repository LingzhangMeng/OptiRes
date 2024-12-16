#' Calculate Average Silhouette Scores Across Resolutions
#'
#' This function computes the average silhouette scores for clustering results across a range of resolutions using a Seurat object. It aids in determining the optimal clustering resolution for single-cell RNA-seq data analysis.
#'
#' @param seu_obj A Seurat object containing preprocessed single-cell data, including UMAP embeddings.
#' @param resolutions Numeric vector of resolution values to evaluate. Default is \code{seq(0.01, 2, 0.01)}.
#'
#' @return A data frame with two columns: \code{Resolution} and \code{SilhouetteScore}, representing the resolution values and their corresponding average silhouette scores.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu_obj' is a preprocessed Seurat object with UMAP embeddings
#' results_df <- scSilhouette(seu_obj, resolutions = seq(0.1, 1.0, 0.1))
#' # Plot the silhouette scores
#' plot(results_df$Resolution, results_df$SilhouetteScore, type = "b",
#'      xlab = "Resolution", ylab = "Average Silhouette Score",
#'      main = "Silhouette Scores Across Resolutions")
#' }
scSilhouette <- function(seu_obj, resolutions = seq(0.01, 2, 0.01)) {
  # Check if UMAP embeddings exist
  if (!"umap" %in% Reductions(seu_obj)) {
    stop("UMAP embeddings not found in the Seurat object. Please run RunUMAP first.")
  }

  silhouette_scores <- numeric(length(resolutions))

  for (i in seq_along(resolutions)) {
    resolution <- resolutions[i]
    # Create a copy to avoid modifying the original object
    seu_tmp <- seu_obj
    # Perform clustering
    seu_tmp <- FindClusters(seu_tmp, resolution = resolution, verbose = FALSE)
    # Get the UMAP coordinates
    umap_coords <- Embeddings(seu_tmp, "umap")
    # Compute silhouette scores
    sil <- silhouette(as.numeric(seu_tmp$seurat_clusters), dist(umap_coords))
    silhouette_scores[i] <- mean(sil[, 3])  # Average silhouette score
  }

  # Create a data frame for plotting
  results_df <- data.frame(Resolution = resolutions, SilhouetteScore = silhouette_scores)

  highest_silhouette_scores <- max(results_df$SilhouetteScore)

  cat("Highest Silhouette Score:", highest_silhouette_scores, "\n")

  # Identify the optimal resolution
  optimal_resolution <- resolutions[which.max(silhouette_scores)]

  cat("Optimal Resolution:", optimal_resolution, "\n")
  return(results_df)
}
