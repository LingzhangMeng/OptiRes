#' Calculate Average Silhouette Scores Across Resolutions
#'
#' This function computes the average silhouette scores for clustering results across a range of resolutions using a Seurat object. It aids in determining the optimal clustering resolution for single-cell RNA-seq data analysis.
#'
#' @param seu_obj A Seurat object containing preprocessed single-cell data, including UMAP embeddings.
#' @param resolutions Numeric vector of resolution values to evaluate. Must be positive. Default is \code{seq(0.01, 2, 0.01)}.
#' @param verbose Logical. Whether to print progress and results. Default is \code{TRUE}.
#'
#' @return A data frame with two columns: \code{Resolution} and \code{SilhouetteScore}, representing the resolution values and their corresponding average silhouette scores.
#' @export
#'
#' @importFrom Seurat FindClusters Embeddings Reductions
#' @importFrom cluster silhouette
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu_obj' is a preprocessed Seurat object with UMAP embeddings
#' results_df <- scSilhouette(seu_obj, resolutions = seq(0.1, 1.0, 0.1))
#' # Plot the silhouette scores
#' Plot_res(results_df)
#' }
scSilhouette <- function(seu_obj, resolutions = seq(0.01, 2, 0.01), verbose = TRUE) {
  # Input validation
  if (!inherits(seu_obj, "Seurat")) {
    stop("The 'seu_obj' parameter must be a valid Seurat object.")
  }
  if (!is.numeric(resolutions) || any(resolutions <= 0)) {
    stop("The 'resolutions' parameter must be a numeric vector with positive values.")
  }
  if (!is.logical(verbose)) {
    stop("The 'verbose' parameter must be logical (TRUE or FALSE).")
  }
  if (!"umap" %in% Reductions(seu_obj)) {
    stop("UMAP embeddings not found in the Seurat object. Please run RunUMAP first.")
  }

  # Initialize results vector
  silhouette_scores <- numeric(length(resolutions))

  if (verbose) cat("Calculating silhouette scores for", length(resolutions), "resolutions...\n")

  # Compute silhouette scores for each resolution
  for (i in seq_along(resolutions)) {
    resolution <- resolutions[i]

    if (verbose) cat(sprintf("Computing resolution at %.2f\n", resolution))

    # Perform clustering without modifying the original object
    seu_tmp <- FindClusters(seu_obj, resolution = resolution, verbose = FALSE)

    # Get UMAP coordinates
    umap_coords <- Embeddings(seu_tmp, "umap")

    # Ensure clusters are valid for silhouette calculation
    clusters <- as.numeric(seu_tmp$seurat_clusters)
    if (length(unique(clusters)) < 2) {
      if (verbose) cat("  Warning: Only one cluster found at resolution", resolution, ". Assigning silhouette score of 0.\n")
      silhouette_scores[i] <- 0
      next
    }

    # Compute silhouette scores
    sil <- silhouette(clusters, dist(umap_coords))
    silhouette_scores[i] <- mean(sil[, 3], na.rm = TRUE)

    if (verbose) cat(sprintf("Silhouette score = %.4f\n", silhouette_scores[i]))
  }

  # Create results data frame
  results_df <- data.frame(Resolution = resolutions, SilhouetteScore = silhouette_scores)

  # Identify the highest silhouette score and optimal resolution
  highest_silhouette_score <- max(results_df$SilhouetteScore, na.rm = TRUE)
  optimal_resolution <- results_df$Resolution[which.max(results_df$SilhouetteScore)]

  if (verbose) {
    cat(sprintf("Highest Silhouette Score: %.4f\n", highest_silhouette_score))
    cat(sprintf("Optimal Resolution: %.2f\n", optimal_resolution))
  }

  return(results_df)
}


#' Plot Silhouette Scores Across Resolutions
#'
#' This function visualizes the silhouette scores for different clustering resolutions and highlights the optimal resolution using ggplot2.
#'
#' @param results_df A data frame with columns \code{Resolution} and \code{SilhouetteScore}, typically the output of \code{scSilhouette}.
#'
#' @return A ggplot object displaying the silhouette scores with the optimal resolution highlighted.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_text labs theme_minimal theme element_blank element_rect element_text
#'
#' @examples
#' \dontrun{
#' # Assuming 'results_df' is the output of scSilhouette
#' p <- Plot_res(results_df)
#' print(p)
#' }
Plot_res <- function(results_df) {
  # Input validation
  if (!is.data.frame(results_df) || !all(c("Resolution", "SilhouetteScore") %in% colnames(results_df))) {
    stop("The 'results_df' parameter must be a data frame with 'Resolution' and 'SilhouetteScore' columns.")
  }
  if (nrow(results_df) == 0) {
    stop("The 'results_df' parameter must not be empty.")
  }

  # Identify the optimal resolution
  optimal_idx <- which.max(results_df$SilhouetteScore)
  optimal_resolution <- results_df$Resolution[optimal_idx]
  optimal_score <- results_df$SilhouetteScore[optimal_idx]

  # Create the plot
  p <- ggplot(results_df, aes(x = Resolution, y = SilhouetteScore)) +
    geom_line(color = "black", linewidth = 0.5) +
    geom_point(color = "black", size = 1.5) +
    geom_point(
      data = results_df[optimal_idx, , drop = FALSE],
      aes(x = Resolution, y = SilhouetteScore),
      color = "red", size = 3
    ) +
    geom_text(
      data = results_df[optimal_idx, , drop = FALSE],
      aes(label = sprintf("Optimal: %.2f", Resolution)),
      vjust = -1, color = "red", size = 4
    ) +
    labs(
      title = "Silhouette Scores Across Resolutions",
      x = "Resolution",
      y = "Average Silhouette Score"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  return(p)
}


#' Plot Colorful Cluster Tree
#'
#' This function generates a colorful dendrogram for the Seurat object's cluster tree, highlighting clusters at the optimal resolution.
#'
#' @param seu_obj A Seurat object containing preprocessed single-cell data.
#' @param results_df A data frame containing the silhouette scores and resolutions, typically the output of \code{scSilhouette}.
#' @param dims Numeric vector specifying the dimensions to use for building the cluster tree. Default is \code{1:5}.
#'
#' @return A dendrogram object with branches and labels colored by cluster.
#' @export
#'
#' @importFrom Seurat FindClusters BuildClusterTree
#' @importFrom stats as.dendrogram
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu_obj' is a preprocessed Seurat object and 'results_df' is the output of scSilhouette
#' results_df <- scSilhouette(seu_obj)
#' dendro <- Plot_ColorfulClusterTree(seu_obj, results_df, dims = 1:10)
#' plot(dendro)
#' }
Plot_ColorfulClusterTree <- function(seu_obj, results_df, dims = 1:5) {
  # Input validation
  if (!inherits(seu_obj, "Seurat")) {
    stop("The 'seu_obj' parameter must be a valid Seurat object.")
  }
  if (!is.data.frame(results_df) || !all(c("Resolution", "SilhouetteScore") %in% colnames(results_df))) {
    stop("The 'results_df' parameter must be a data frame with 'Resolution' and 'SilhouetteScore' columns.")
  }
  if (!is.numeric(dims) || any(dims < 1)) {
    stop("The 'dims' parameter must be a numeric vector with positive integers.")
  }
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("The 'dendextend' package is required but not installed. Please install it using install.packages('dendextend').")
  }

  # Identify the optimal resolution
  optimal_resolution <- results_df$Resolution[which.max(results_df$SilhouetteScore)]

  # Perform clustering with the optimal resolution
  seu_obj <- FindClusters(seu_obj, resolution = optimal_resolution, verbose = FALSE)

  # Build the cluster tree
  seu_obj <- BuildClusterTree(seu_obj, dims = dims)

  # Check if cluster tree was successfully built
  if (is.null(seu_obj@tools$BuildClusterTree)) {
    stop("Failed to build cluster tree. Ensure sufficient data and valid dimensions.")
  }

  # Convert to dendrogram
  dendro <- as.dendrogram(seu_obj@tools$BuildClusterTree)

  # Get number of clusters
  num_clusters <- length(unique(seu_obj@active.ident))

  # Color the dendrogram
  dendro_colored <- dendextend::set(dendro, "branches_k_color", k = num_clusters) %>%
    dendextend::set("labels_colors", k = num_clusters) %>%
    dendextend::set("labels_cex", 0.8)

  # Plot the dendrogram
  plot(dendro_colored, main = sprintf("Cluster Tree at Resolution %.2f", optimal_resolution))

  return(dendro_colored)
}
