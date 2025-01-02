#' Calculate Average Silhouette Scores Across Resolutions
#'
#' This function computes the average silhouette scores for clustering results across a range of resolutions using a Seurat object. It aids in determining the optimal clustering resolution for single-cell RNA-seq data analysis.
#'
#' @param seu_obj A Seurat object containing preprocessed single-cell data, including UMAP embeddings.
#' @param resolutions Numeric vector of resolution values to evaluate. Default is \code{seq(0.01, 2, 0.01)}.
#' @param verbose Logical. Whether to print progress and results. Default is \code{TRUE}.
#'
#' @return A data frame with two columns: \code{Resolution} and \code{SilhouetteScore}, representing the resolution values and their corresponding average silhouette scores.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu_obj' is a preprocessed Seurat object with UMAP embeddings
#' results_df <- scSilhouette(seu_obj, resolutions = seq(0.1, 1.0, 0.1))
#' # Plot the silhouette scores
#' results_df <- scSilhouette(seu_obj)
#' Plot_res(results_df)
#' }
scSilhouette <- function(seu_obj, resolutions = seq(0.01, 2, 0.01), verbose = TRUE) {
  # Check if UMAP embeddings exist
  if (!"umap" %in% Reductions(seu_obj)) {
    stop("UMAP embeddings not found in the Seurat object. Please run RunUMAP first.")
  }

  silhouette_scores <- numeric(length(resolutions))

  if (verbose) cat("Calculating silhouette scores for", length(resolutions), "resolutions...\n")

  for (i in seq_along(resolutions)) {
    resolution <- resolutions[i]

    if (verbose) cat("Calculating silhouette score at resolution:", resolution, "\n")

    # Create a copy to avoid modifying the original object
    seu_tmp <- seu_obj

    # Perform clustering
    seu_tmp <- FindClusters(seu_tmp, resolution = resolution, verbose = FALSE)

    # Get the UMAP coordinates
    umap_coords <- Embeddings(seu_tmp, "umap")

    # Compute silhouette scores
    sil <- silhouette(as.numeric(seu_tmp$seurat_clusters), dist(umap_coords))
    silhouette_scores[i] <- mean(sil[, 3])  # Average silhouette score

    if (verbose) cat("  Silhouette score =", silhouette_scores[i], "\n")
  }

  # Create a data frame for the results
  results_df <- data.frame(Resolution = resolutions, SilhouetteScore = silhouette_scores)

  # Identify the highest silhouette score and the optimal resolution
  highest_silhouette_score <- max(results_df$SilhouetteScore)
  optimal_resolution <- resolutions[which.max(silhouette_scores)]

  if (verbose) {
    cat("Highest Silhouette Score:", highest_silhouette_score, "\n")
    cat("Optimal Resolution:", optimal_resolution, "\n")
  }

  return(results_df)
}


#' Plot Silhouette Scores Across Resolutions
#'
#' This function visualizes the silhouette scores for different clustering resolutions and highlights the optimal resolution.
#'
#' @param results_df A data frame with columns \code{Resolution} and \code{SilhouetteScore}, typically the output of \code{scSilhouette}.
#'
#' @return A ggplot object showing the silhouette scores and the optimal resolution.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'results_df' is the output of scSilhouette
#' Plot_res(results_df)
#' }
Plot_res <- function(results_df) {
  # Validate input
  if (!all(c("Resolution", "SilhouetteScore") %in% colnames(results_df))) {
    stop("The input data frame must contain 'Resolution' and 'SilhouetteScore' columns.")
  }

  # Identify the optimal resolution
  optimal_resolution <- results_df$Resolution[which.max(results_df$SilhouetteScore)]
  optimal_score <- max(results_df$SilhouetteScore)

  # Create the plot
  p <- ggplot(results_df, aes(x = Resolution, y = SilhouetteScore)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_point(data = results_df[results_df$Resolution == optimal_resolution, ],
               aes(x = Resolution, y = SilhouetteScore),
               color = "red", size = 3) +  # Highlight optimal resolution
    geom_text(data = results_df[results_df$Resolution == optimal_resolution, ],
              aes(label = paste("Optimal:", round(Resolution, 2))),
              vjust = -1, color = "red", size = 4) +  # Add label for optimal resolution
    labs(title = "Silhouette Scores for Different Resolutions",
         x = "Resolution",
         y = "Average Silhouette Score") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA),
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}


#' Plot Colorful Cluster Tree
#'
#' This function generates a colorful dendrogram for the Seurat object's cluster tree.
#'
#' @param seu_obj A Seurat object containing preprocessed single-cell data.
#' @param results_df A data frame containing the silhouette scores and resolutions, typically the output of \code{scSilhouette}.
#'
#' @return A dendrogram object with branches and labels colored by cluster.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu_obj' is a preprocessed Seurat object and 'results_df' is the output of scSilhouette
#' results_df <- scSilhouette(seu_obj)
#' Plot_ColorfulClusterTree(seu_obj, results_df)
#' }
Plot_ColorfulClusterTree <- function(seu_obj, results_df) {
  # Validate input
  if (!all(c("Resolution", "SilhouetteScore") %in% colnames(results_df))) {
    stop("The input results_df must contain 'Resolution' and 'SilhouetteScore' columns.")
  }

  # Identify the optimal resolution
  optimal_resolution <- results_df$Resolution[which.max(results_df$SilhouetteScore)]

  # Perform clustering with the optimal resolution
  seu_obj <- FindClusters(seu_obj, resolution = optimal_resolution)

  # Build the cluster tree
  seu_obj <- BuildClusterTree(seu_obj, dims = 1:5)  # Adjust `dims` as needed

  # Convert the cluster tree to a dendrogram
  dendro <- as.dendrogram(seu_obj@tools$BuildClusterTree)

  # Set the number of clusters
  num_clusters <- length(unique(seu_obj@active.ident))

  # Color the dendrogram branches and labels
  library(dendextend)  # Ensure the dendextend package is installed
  dendro_colored <- dendro %>%
    set("branches_k_color", k = num_clusters) %>%
    set("labels_colors", k = num_clusters) %>%
    set("labels_cex", 0.8)

  # Plot the colorful dendrogram
  plot(dendro_colored, main = "Cluster Tree with Colored Branches and Labels")

  # Return the dendrogram object (optional)
  return(dendro_colored)
}
