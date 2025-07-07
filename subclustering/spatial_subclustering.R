#!/usr/bin/env Rscript

# ============================================================================
# Spatial-Aware Subclustering Analysis
# ============================================================================
# 
# This script performs subclustering analysis that incorporates spatial
# information to identify spatially coherent cell populations and validate
# clustering results using spatial patterns.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(Matrix)
  library(sp)
  library(spdep)
})

# Parameters
CLUSTERING_RESOLUTIONS <- seq(0.1, 2.0, 0.1)  # Clustering resolutions to test
PCA_DIMS <- 1:30                              # PCA dimensions to use
SPATIAL_DIMS <- 1:2                           # Spatial dimensions to use
SPATIAL_WEIGHT <- 0.3                         # Weight for spatial information in clustering
MIN_SPATIAL_COHERENCE <- 0.5                  # Minimum spatial coherence score

# Spatial analysis parameters
SPATIAL_RADIUS <- 100                         # Radius for spatial neighborhood analysis
MIN_NEIGHBORS <- 5                            # Minimum neighbors for spatial analysis
MORAN_SIGNIFICANCE <- 0.05                    # P-value threshold for Moran's I test

# Input/Output paths
INPUT_SEURAT <- "../data/spatial_integrated_object.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Helper Functions
# ============================================================================

# Function to calculate spatial coherence score
calculate_spatial_coherence <- function(coords, clusters, k = 10) {
  if (nrow(coords) < k) {
    return(0)
  }
  
  # Find k nearest neighbors for each cell
  knn_indices <- RANN::nn2(coords, k = k + 1)$nn.idx[, -1]  # Exclude self
  
  # Calculate coherence for each cell
  coherence_scores <- sapply(1:nrow(coords), function(i) {
    neighbors <- knn_indices[i, ]
    same_cluster <- clusters[neighbors] == clusters[i]
    mean(same_cluster)
  })
  
  return(mean(coherence_scores))
}

# Function to calculate Moran's I for spatial autocorrelation
calculate_morans_i <- function(coords, values, radius = 100) {
  if (length(unique(values)) == 1) {
    return(list(morans_i = 0, p_value = 1))
  }
  
  # Create spatial weights matrix
  tryCatch({
    # Convert to spatial points
    coords_sp <- sp::SpatialPoints(coords)
    
    # Create distance-based neighbors
    nb <- spdep::dnearneigh(coords_sp, 0, radius)
    
    # Create weights
    if (length(nb[[1]]) == 0) {
      return(list(morans_i = 0, p_value = 1))
    }
    
    weights <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    
    # Calculate Moran's I
    moran_result <- spdep::moran.test(values, weights, zero.policy = TRUE)
    
    return(list(
      morans_i = moran_result$estimate["Moran I statistic"],
      p_value = moran_result$p.value
    ))
    
  }, error = function(e) {
    return(list(morans_i = 0, p_value = 1))
  })
}

# Function to identify spatially compact clusters
identify_compact_clusters <- function(coords, clusters, threshold = 0.5) {
  cluster_coherence <- sapply(unique(clusters), function(cl) {
    cluster_coords <- coords[clusters == cl, , drop = FALSE]
    if (nrow(cluster_coords) < 3) return(0)
    
    # Calculate spatial compactness using convex hull area
    if (nrow(cluster_coords) >= 3) {
      tryCatch({
        hull <- chull(cluster_coords)
        hull_area <- sp::Polygon(cluster_coords[hull, ])@area
        n_cells <- nrow(cluster_coords)
        compactness <- n_cells / (hull_area + 1)  # Add 1 to avoid division by zero
        return(compactness)
      }, error = function(e) {
        return(0)
      })
    } else {
      return(0)
    }
  })
  
  names(cluster_coherence) <- unique(clusters)
  return(cluster_coherence)
}

# Function to create spatial clustering plot
create_spatial_plot <- function(coords, clusters, title = "Spatial Clusters", 
                               point_size = 0.5, alpha = 0.8) {
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    cluster = as.factor(clusters)
  )
  
  # Use a color palette that works for many clusters
  n_clusters <- length(unique(clusters))
  if (n_clusters <= 12) {
    colors <- RColorBrewer::brewer.pal(min(n_clusters, 12), "Set3")
  } else {
    colors <- rainbow(n_clusters)
  }
  
  ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
    geom_point(size = point_size, alpha = alpha) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank()
    ) +
    ggtitle(title) +
    scale_color_manual(values = colors)
}

# Function to perform spatial-aware clustering
spatial_aware_clustering <- function(sobj, spatial_coords, spatial_weight = 0.3, 
                                   resolution = 0.5, dims = 1:30) {
  
  # Standard clustering
  sobj <- FindNeighbors(sobj, dims = dims, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = resolution, verbose = FALSE)
  standard_clusters <- Idents(sobj)
  
  # Calculate spatial weights
  spatial_dist <- as.matrix(dist(spatial_coords))
  spatial_sim <- exp(-spatial_dist / median(spatial_dist))
  
  # Get expression similarity from SNN graph
  if ("RNA_snn" %in% names(sobj@graphs)) {
    expr_sim <- as.matrix(sobj@graphs$RNA_snn)
  } else {
    # Fallback: use correlation
    expr_data <- GetAssayData(sobj, slot = "scale.data")
    expr_sim <- cor(t(expr_data))
  }
  
  # Combine expression and spatial similarity
  combined_sim <- (1 - spatial_weight) * expr_sim + spatial_weight * spatial_sim
  
  # Perform clustering on combined similarity
  combined_dist <- as.dist(1 - combined_sim)
  spatial_clusters <- cutree(hclust(combined_dist), k = length(unique(standard_clusters)))
  
  return(list(
    standard_clusters = standard_clusters,
    spatial_clusters = spatial_clusters,
    combined_similarity = combined_sim
  ))
}

# ============================================================================
# Load Data and Prepare Spatial Information
# ============================================================================

cat("Loading spatial transcriptomics data...\\n")

if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Extract spatial coordinates
spatial_coords <- NULL
coord_source <- "unknown"

# Try different spatial coordinate sources
if ("Spatial" %in% names(sobj@reductions)) {
  spatial_coords <- Embeddings(sobj, "Spatial")
  coord_source <- "Spatial reduction"
} else if ("SpatialCenter" %in% names(sobj@reductions)) {
  spatial_coords <- Embeddings(sobj, "SpatialCenter")
  coord_source <- "SpatialCenter reduction"
} else if (all(c("x", "y") %in% colnames(sobj@meta.data))) {
  spatial_coords <- as.matrix(sobj@meta.data[, c("x", "y")])
  coord_source <- "metadata coordinates"
} else if (all(c("imagerow", "imagecol") %in% colnames(sobj@meta.data))) {
  spatial_coords <- as.matrix(sobj@meta.data[, c("imagecol", "imagerow")])
  coord_source <- "image coordinates"
} else {
  stop("No spatial coordinates found in the Seurat object")
}

cat("Using spatial coordinates from:", coord_source, "\\n")
cat("Coordinate range: X [", round(min(spatial_coords[,1]), 2), ", ", 
    round(max(spatial_coords[,1]), 2), "], Y [", round(min(spatial_coords[,2]), 2), 
    ", ", round(max(spatial_coords[,2]), 2), "]\\n")

# Check for multiple samples/slices
sample_col <- "sample_id"
if (!"sample_id" %in% colnames(sobj@meta.data)) {
  if ("id" %in% colnames(sobj@meta.data)) {
    sample_col <- "id"
  } else if ("orig.ident" %in% colnames(sobj@meta.data)) {
    sample_col <- "orig.ident"
  }
}

n_samples <- length(unique(sobj@meta.data[[sample_col]]))
cat("Number of samples/slices:", n_samples, "\\n")

# ============================================================================
# Standard Multi-Resolution Clustering
# ============================================================================

cat("Performing standard multi-resolution clustering...\\n")

# Ensure PCA and neighbors are computed
if (!"pca" %in% names(sobj@reductions)) {
  sobj <- RunPCA(sobj, verbose = FALSE)
}

if (!"RNA_snn" %in% names(sobj@graphs)) {
  sobj <- FindNeighbors(sobj, dims = PCA_DIMS, verbose = FALSE)
}

# Perform clustering at multiple resolutions
sobj <- FindClusters(sobj, resolution = CLUSTERING_RESOLUTIONS, verbose = FALSE)

# Extract clustering results
resolution_columns <- paste0("RNA_snn_res.", CLUSTERING_RESOLUTIONS)
available_columns <- resolution_columns[resolution_columns %in% colnames(sobj@meta.data)]

cat("Generated", length(available_columns), "clustering results\\n")

# ============================================================================
# Spatial Clustering Analysis
# ============================================================================

cat("Performing spatial-aware clustering analysis...\\n")

# Initialize results
spatial_results <- data.frame(
  resolution = numeric(),
  n_clusters = numeric(),
  spatial_coherence = numeric(),
  morans_i = numeric(),
  morans_p = numeric(),
  compact_clusters = numeric(),
  stringsAsFactors = FALSE
)

# Test spatial-aware clustering at different resolutions
for (i in seq_along(CLUSTERING_RESOLUTIONS)) {
  res <- CLUSTERING_RESOLUTIONS[i]
  res_col <- paste0("RNA_snn_res.", res)
  
  if (!res_col %in% colnames(sobj@meta.data)) next
  
  clusters <- sobj@meta.data[[res_col]]
  n_clusters <- length(unique(clusters))
  
  cat("Analyzing resolution", res, "(", n_clusters, "clusters)...\\n")
  
  # Calculate spatial coherence
  coherence <- calculate_spatial_coherence(spatial_coords, clusters)
  
  # Calculate Moran's I for spatial autocorrelation
  moran_result <- calculate_morans_i(spatial_coords, as.numeric(as.factor(clusters)), 
                                    radius = SPATIAL_RADIUS)
  
  # Identify spatially compact clusters
  compact_scores <- identify_compact_clusters(spatial_coords, clusters)
  n_compact <- sum(compact_scores > quantile(compact_scores, 0.5, na.rm = TRUE), na.rm = TRUE)
  
  # Store results
  spatial_results <- rbind(spatial_results, data.frame(
    resolution = res,
    n_clusters = n_clusters,
    spatial_coherence = coherence,
    morans_i = moran_result$morans_i,
    morans_p = moran_result$p_value,
    compact_clusters = n_compact,
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# Spatial-Aware Clustering with Optimal Parameters
# ============================================================================

cat("Performing optimized spatial-aware clustering...\\n")

# Find resolution with best spatial coherence
optimal_idx <- which.max(spatial_results$spatial_coherence)
optimal_resolution <- spatial_results$resolution[optimal_idx]

cat("Optimal spatial resolution:", optimal_resolution, "\\n")

# Perform spatial-aware clustering at optimal resolution
spatial_clustering_result <- spatial_aware_clustering(
  sobj, spatial_coords, 
  spatial_weight = SPATIAL_WEIGHT, 
  resolution = optimal_resolution, 
  dims = PCA_DIMS
)

# Add spatial clustering results to object
sobj$spatial_clusters <- spatial_clustering_result$spatial_clusters
sobj$standard_clusters_optimal <- spatial_clustering_result$standard_clusters

# ============================================================================
# Visualization
# ============================================================================

cat("Creating visualizations...\\n")

# 1. Spatial coherence vs resolution
coherence_plot <- ggplot(spatial_results, aes(x = resolution, y = spatial_coherence)) +
  geom_line(size = 1, color = "blue") +
  geom_point(size = 2, color = "blue") +
  geom_vline(xintercept = optimal_resolution, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Spatial Coherence vs Resolution",
       x = "Resolution", y = "Spatial Coherence Score") +
  annotate("text", x = optimal_resolution + 0.1, y = max(spatial_results$spatial_coherence) * 0.9,
           label = paste("Optimal:", optimal_resolution), color = "red")

# 2. Number of clusters vs spatial metrics
metrics_plot <- spatial_results %>%
  select(resolution, n_clusters, spatial_coherence, compact_clusters) %>%
  pivot_longer(cols = c(n_clusters, compact_clusters), 
               names_to = "metric", values_to = "count") %>%
  ggplot(aes(x = resolution, y = count, color = metric)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Clustering Metrics vs Resolution",
       x = "Resolution", y = "Count") +
  scale_color_discrete(name = "Metric", 
                      labels = c("Compact Clusters", "Total Clusters"))

# 3. Spatial plots comparison
standard_spatial_plot <- create_spatial_plot(
  spatial_coords, 
  spatial_clustering_result$standard_clusters,
  title = "Standard Clustering"
)

spatial_aware_plot <- create_spatial_plot(
  spatial_coords, 
  spatial_clustering_result$spatial_clusters,
  title = "Spatial-Aware Clustering"
)

# 4. UMAP comparison
if ("umap" %in% names(sobj@reductions)) {
  umap_standard <- DimPlot(sobj, group.by = "standard_clusters_optimal", 
                          label = TRUE, pt.size = 0.5) +
    ggtitle("Standard Clustering (UMAP)")
  
  umap_spatial <- DimPlot(sobj, group.by = "spatial_clusters", 
                         label = TRUE, pt.size = 0.5) +
    ggtitle("Spatial-Aware Clustering (UMAP)")
} else {
  # Run UMAP if not available
  sobj <- RunUMAP(sobj, dims = PCA_DIMS, verbose = FALSE)
  
  umap_standard <- DimPlot(sobj, group.by = "standard_clusters_optimal", 
                          label = TRUE, pt.size = 0.5) +
    ggtitle("Standard Clustering (UMAP)")
  
  umap_spatial <- DimPlot(sobj, group.by = "spatial_clusters", 
                         label = TRUE, pt.size = 0.5) +
    ggtitle("Spatial-Aware Clustering (UMAP)")
}

# Save plots
ggsave(file.path(PLOTS_DIR, "spatial_coherence_analysis.pdf"), 
       (coherence_plot | metrics_plot), width = 16, height = 6)

ggsave(file.path(PLOTS_DIR, "spatial_clustering_comparison.pdf"), 
       (standard_spatial_plot | spatial_aware_plot), width = 16, height = 8)

ggsave(file.path(PLOTS_DIR, "umap_clustering_comparison.pdf"), 
       (umap_standard | umap_spatial), width = 16, height = 8)

# ============================================================================
# Cluster Validation and Characterization
# ============================================================================

cat("Validating and characterizing spatial clusters...\\n")

# Compare standard vs spatial-aware clustering
cluster_comparison <- data.frame(
  cell_id = colnames(sobj),
  standard_cluster = spatial_clustering_result$standard_clusters,
  spatial_cluster = spatial_clustering_result$spatial_clusters,
  x = spatial_coords[, 1],
  y = spatial_coords[, 2]
)

# Calculate agreement between methods
agreement <- mean(cluster_comparison$standard_cluster == cluster_comparison$spatial_cluster)
cat("Clustering agreement:", round(100 * agreement, 1), "%\\n")

# Find markers for spatial clusters
Idents(sobj) <- "spatial_clusters"
spatial_markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, 
                                 logfc.threshold = 0.25, verbose = FALSE)

# Cluster summary with spatial metrics
cluster_summary <- cluster_comparison %>%
  group_by(spatial_cluster) %>%
  summarise(
    n_cells = n(),
    x_mean = mean(x),
    y_mean = mean(y),
    x_sd = sd(x),
    y_sd = sd(y),
    spatial_spread = sqrt(var(x) + var(y)),
    .groups = 'drop'
  )

# Add spatial coherence for each cluster
cluster_summary$coherence <- sapply(cluster_summary$spatial_cluster, function(cl) {
  cluster_cells <- cluster_comparison$spatial_cluster == cl
  cluster_coords <- spatial_coords[cluster_cells, , drop = FALSE]
  cluster_labels <- rep(1, sum(cluster_cells))
  
  if (sum(cluster_cells) < 3) return(1)
  
  calculate_spatial_coherence(cluster_coords, cluster_labels, k = min(5, sum(cluster_cells) - 1))
})

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating spatial subclustering report...\\n")

# Summary statistics
n_cells <- ncol(sobj)
n_standard_clusters <- length(unique(spatial_clustering_result$standard_clusters))
n_spatial_clusters <- length(unique(spatial_clustering_result$spatial_clusters))
optimal_coherence <- spatial_results$spatial_coherence[optimal_idx]

report_text <- paste0(
  "Spatial-Aware Subclustering Analysis Report\\n",
  "==========================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Dataset: ", n_cells, " cells from ", n_samples, " samples\\n",
  "Spatial coordinates: ", coord_source, "\\n\\n",
  "Analysis Parameters:\\n",
  "  - Resolution range: ", min(CLUSTERING_RESOLUTIONS), " - ", max(CLUSTERING_RESOLUTIONS), "\\n",
  "  - PCA dimensions: ", paste(range(PCA_DIMS), collapse = "-"), "\\n",
  "  - Spatial weight: ", SPATIAL_WEIGHT, "\\n",
  "  - Spatial radius: ", SPATIAL_RADIUS, "\\n\\n",
  "Results Summary:\\n",
  "  - Optimal resolution: ", optimal_resolution, "\\n",
  "  - Standard clusters: ", n_standard_clusters, "\\n",
  "  - Spatial-aware clusters: ", n_spatial_clusters, "\\n",
  "  - Clustering agreement: ", round(100 * agreement, 1), "%\\n",
  "  - Optimal spatial coherence: ", round(optimal_coherence, 3), "\\n\\n",
  "Spatial Validation:\\n",
  "  - Moran's I: ", round(spatial_results$morans_i[optimal_idx], 3), "\\n",
  "  - Moran's I p-value: ", format(spatial_results$morans_p[optimal_idx], scientific = TRUE), "\\n",
  "  - Spatially compact clusters: ", spatial_results$compact_clusters[optimal_idx], "\\n\\n"
)

# Add cluster-specific information
report_text <- paste0(report_text, "Cluster Characteristics:\\n")
for (i in 1:nrow(cluster_summary)) {
  cluster_info <- cluster_summary[i, ]
  report_text <- paste0(report_text,
    "  Cluster ", cluster_info$spatial_cluster, ":\\n",
    "    - Cells: ", cluster_info$n_cells, "\\n",
    "    - Spatial coherence: ", round(cluster_info$coherence, 3), "\\n",
    "    - Spatial spread: ", round(cluster_info$spatial_spread, 1), "\\n\\n"
  )
}

report_text <- paste0(report_text,
  "Output Files:\\n",
  "  - Spatial analysis results: spatial_clustering_results.csv\\n",
  "  - Cluster summary: spatial_cluster_summary.csv\\n",
  "  - Marker genes: spatial_cluster_markers.csv\\n",
  "  - Processed object: spatial_subclustered_object.rds\\n\\n",
  "Visualization Files:\\n",
  "  - spatial_coherence_analysis.pdf: Resolution analysis\\n",
  "  - spatial_clustering_comparison.pdf: Spatial comparison\\n",
  "  - umap_clustering_comparison.pdf: UMAP comparison\\n\\n",
  "Recommendations:\\n",
  "1. Use spatial-aware clustering for spatially organized tissues\\n",
  "2. Validate clusters with known spatial organization patterns\\n",
  "3. Consider higher spatial weights for more spatial influence\\n",
  "4. Review clusters with low spatial coherence\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "spatial_subclustering_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\\n")

# Save spatial analysis results
write.csv(spatial_results, file.path(OUTPUT_DIR, "spatial_clustering_results.csv"), 
          row.names = FALSE)

# Save cluster summary
write.csv(cluster_summary, file.path(OUTPUT_DIR, "spatial_cluster_summary.csv"), 
          row.names = FALSE)

# Save cluster comparison
write.csv(cluster_comparison, file.path(OUTPUT_DIR, "clustering_method_comparison.csv"), 
          row.names = FALSE)

# Save marker genes
if (nrow(spatial_markers) > 0) {
  write.csv(spatial_markers, file.path(OUTPUT_DIR, "spatial_cluster_markers.csv"), 
            row.names = FALSE)
}

# Save processed Seurat object
saveRDS(sobj, file.path(OUTPUT_DIR, "spatial_subclustered_object.rds"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_spatial_subclustering.txt"))

cat("\\nSpatial-aware subclustering analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nKey findings:\\n")
cat("- Optimal resolution:", optimal_resolution, "\\n")
cat("- Spatial coherence score:", round(optimal_coherence, 3), "\\n")
cat("- Clustering method agreement:", round(100 * agreement, 1), "%\\n")
cat("- Spatial clusters identified:", n_spatial_clusters, "\\n")

if (nrow(spatial_markers) > 0) {
  cat("- Spatial marker genes:", nrow(spatial_markers), "\\n")
}

cat("\\nNext steps:\\n")
cat("1. Review spatial clustering patterns\\n")
cat("2. Validate with known tissue organization\\n")
cat("3. Use for spatially-informed downstream analysis\\n")