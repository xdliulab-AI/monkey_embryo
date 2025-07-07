#!/usr/bin/env Rscript

# ============================================================================
# Multi-Resolution Clustering Analysis
# ============================================================================
# 
# This script performs comprehensive multi-resolution clustering analysis
# to identify optimal clustering parameters and explore hierarchical
# relationships within cell populations.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(clustree)
  library(RColorBrewer)
  library(pheatmap)
  library(Matrix)
})

# Parameters
RESOLUTION_RANGE <- seq(0.1, 3.0, 0.1)  # Range of resolutions to test
PCA_DIMS <- 1:30                        # PCA dimensions to use
METRICS_TO_CALCULATE <- c("silhouette", "stability", "modularity")  # Clustering metrics
MIN_CLUSTER_SIZE <- 10                  # Minimum cells per cluster for analysis

# Input/Output paths
INPUT_SEURAT <- "../data/processed_seurat_object.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Helper Functions
# ============================================================================

# Function to calculate silhouette scores
calculate_silhouette_score <- function(sobj, reduction = "pca", dims = 1:30, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- Idents(sobj)
  }
  
  embeddings <- Embeddings(sobj, reduction = reduction)[, dims]
  dist_matrix <- dist(embeddings)
  
  # Calculate silhouette scores
  sil_scores <- cluster::silhouette(as.numeric(as.factor(clusters)), dist_matrix)
  
  return(mean(sil_scores[, 3]))
}

# Function to calculate cluster stability across resolutions
calculate_cluster_stability <- function(cluster_matrix) {
  n_res <- ncol(cluster_matrix)
  stability_scores <- numeric(n_res - 1)
  
  for (i in 1:(n_res - 1)) {
    # Calculate adjusted rand index between consecutive resolutions
    ari <- mclust::adjustedRandIndex(cluster_matrix[, i], cluster_matrix[, i + 1])
    stability_scores[i] <- ari
  }
  
  return(mean(stability_scores))
}

# Function to calculate modularity
calculate_modularity <- function(sobj, clusters = NULL) {
  if (is.null(clusters)) {
    clusters <- Idents(sobj)
  }
  
  # Get the SNN graph
  if ("RNA_snn" %in% names(sobj@graphs)) {
    snn_graph <- sobj@graphs$RNA_snn
  } else {
    stop("SNN graph not found. Please run FindNeighbors first.")
  }
  
  # Convert to igraph format
  snn_igraph <- igraph::graph_from_adjacency_matrix(snn_graph, weighted = TRUE, mode = "undirected")
  
  # Calculate modularity
  modularity <- igraph::modularity(snn_igraph, as.numeric(as.factor(clusters)))
  
  return(modularity)
}

# Function to create resolution comparison plot
create_resolution_plot <- function(sobj, resolutions, ncol = 4) {
  plots <- list()
  
  for (i in seq_along(resolutions)) {
    res <- resolutions[i]
    res_col <- paste0("RNA_snn_res.", res)
    
    if (res_col %in% colnames(sobj@meta.data)) {
      n_clusters <- length(unique(sobj@meta.data[[res_col]]))
      
      p <- DimPlot(sobj, group.by = res_col, label = TRUE, pt.size = 0.3) +
        ggtitle(paste0("Resolution ", res, " (", n_clusters, " clusters)")) +
        theme(legend.position = "none",
              plot.title = element_text(size = 10))
      
      plots[[length(plots) + 1]] <- p
    }
  }
  
  return(wrap_plots(plots, ncol = ncol))
}

# Function to analyze cluster composition
analyze_cluster_composition <- function(sobj, metadata_col = "celltype") {
  if (!metadata_col %in% colnames(sobj@meta.data)) {
    return(NULL)
  }
  
  composition_data <- sobj@meta.data %>%
    group_by(seurat_clusters, !!sym(metadata_col)) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(proportion = n_cells / sum(n_cells)) %>%
    ungroup()
  
  return(composition_data)
}

# ============================================================================
# Load Data
# ============================================================================

cat("Loading Seurat object...\\n")

if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Check if clustering has been performed
if (!"RNA_snn" %in% names(sobj@graphs)) {
  cat("SNN graph not found. Computing neighbors...\\n")
  sobj <- FindNeighbors(sobj, dims = PCA_DIMS, verbose = FALSE)
}

# ============================================================================
# Multi-Resolution Clustering
# ============================================================================

cat("Performing multi-resolution clustering...\\n")

# Perform clustering at multiple resolutions
sobj <- FindClusters(sobj, resolution = RESOLUTION_RANGE, verbose = FALSE)

# Extract clustering results
resolution_columns <- paste0("RNA_snn_res.", RESOLUTION_RANGE)
available_columns <- resolution_columns[resolution_columns %in% colnames(sobj@meta.data)]

cat("Generated clustering results for", length(available_columns), "resolutions\\n")

# Create matrix of clustering results
cluster_matrix <- sobj@meta.data[, available_columns]
colnames(cluster_matrix) <- gsub("RNA_snn_res.", "", colnames(cluster_matrix))

# ============================================================================
# Clustering Quality Assessment
# ============================================================================

cat("Assessing clustering quality...\\n")

# Initialize results data frame
quality_metrics <- data.frame(
  resolution = as.numeric(gsub("RNA_snn_res.", "", available_columns)),
  n_clusters = numeric(length(available_columns)),
  n_singletons = numeric(length(available_columns)),
  largest_cluster_prop = numeric(length(available_columns)),
  silhouette_score = numeric(length(available_columns)),
  modularity = numeric(length(available_columns))
)

# Calculate metrics for each resolution
for (i in seq_along(available_columns)) {
  res_col <- available_columns[i]
  clusters <- sobj@meta.data[[res_col]]
  
  # Basic cluster statistics
  cluster_sizes <- table(clusters)
  quality_metrics$n_clusters[i] <- length(cluster_sizes)
  quality_metrics$n_singletons[i] <- sum(cluster_sizes < MIN_CLUSTER_SIZE)
  quality_metrics$largest_cluster_prop[i] <- max(cluster_sizes) / length(clusters)
  
  # Silhouette score
  tryCatch({
    quality_metrics$silhouette_score[i] <- calculate_silhouette_score(sobj, clusters = clusters)
  }, error = function(e) {
    quality_metrics$silhouette_score[i] <- NA
  })
  
  # Modularity
  tryCatch({
    quality_metrics$modularity[i] <- calculate_modularity(sobj, clusters = clusters)
  }, error = function(e) {
    quality_metrics$modularity[i] <- NA
  })
}

# Calculate stability scores
cat("Calculating cluster stability...\\n")
quality_metrics$stability <- NA
for (i in 2:nrow(quality_metrics)) {
  res_prev <- available_columns[i-1]
  res_curr <- available_columns[i]
  
  if (both_exist <- res_prev %in% colnames(sobj@meta.data) && res_curr %in% colnames(sobj@meta.data)) {
    ari <- mclust::adjustedRandIndex(sobj@meta.data[[res_prev]], sobj@meta.data[[res_curr]])
    quality_metrics$stability[i] <- ari
  }
}

# ============================================================================
# Visualization
# ============================================================================

cat("Creating visualizations...\\n")

# 1. Multi-resolution clustering overview
selected_resolutions <- c(0.2, 0.5, 0.8, 1.0, 1.5, 2.0)
selected_resolutions <- selected_resolutions[selected_resolutions %in% quality_metrics$resolution]

if (length(selected_resolutions) > 0) {
  resolution_plot <- create_resolution_plot(sobj, selected_resolutions, ncol = 3)
  ggsave(file.path(PLOTS_DIR, "multi_resolution_overview.pdf"), 
         resolution_plot, width = 15, height = 10)
}

# 2. Quality metrics plots
metrics_plot <- quality_metrics %>%
  pivot_longer(cols = c(silhouette_score, modularity, stability), 
               names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = resolution, y = value, color = metric)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_y", nrow = 3) +
  theme_minimal() +
  labs(title = "Clustering Quality Metrics", 
       x = "Resolution", y = "Score") +
  theme(legend.position = "none")

ggsave(file.path(PLOTS_DIR, "clustering_quality_metrics.pdf"), 
       metrics_plot, width = 10, height = 8)

# 3. Number of clusters vs resolution
cluster_count_plot <- ggplot(quality_metrics, aes(x = resolution, y = n_clusters)) +
  geom_line(size = 1, color = "blue") +
  geom_point(size = 2, color = "blue") +
  geom_line(aes(y = n_singletons), size = 1, color = "red", linetype = "dashed") +
  geom_point(aes(y = n_singletons), size = 2, color = "red") +
  theme_minimal() +
  labs(title = "Cluster Count vs Resolution",
       x = "Resolution", 
       y = "Number of Clusters",
       subtitle = "Blue: Total clusters, Red: Small clusters (<10 cells)") +
  scale_y_continuous(breaks = scales::pretty_breaks())

ggsave(file.path(PLOTS_DIR, "cluster_count_vs_resolution.pdf"), 
       cluster_count_plot, width = 10, height = 6)

# 4. Clustree visualization
if (ncol(cluster_matrix) > 1) {
  tryCatch({
    clustree_plot <- clustree(cluster_matrix, prefix = "", suffix = "") +
      theme(legend.position = "bottom")
    
    ggsave(file.path(PLOTS_DIR, "clustree_resolution_tree.pdf"), 
           clustree_plot, width = 12, height = 10)
  }, error = function(e) {
    cat("Warning: Could not create clustree plot:", e$message, "\\n")
  })
}

# ============================================================================
# Optimal Resolution Selection
# ============================================================================

cat("Identifying optimal resolution...\\n")

# Composite score for resolution selection
# Normalize metrics to 0-1 scale
quality_metrics_norm <- quality_metrics %>%
  mutate(
    silhouette_norm = ifelse(is.na(silhouette_score), 0, 
                            (silhouette_score - min(silhouette_score, na.rm = TRUE)) / 
                            (max(silhouette_score, na.rm = TRUE) - min(silhouette_score, na.rm = TRUE))),
    modularity_norm = ifelse(is.na(modularity), 0,
                            (modularity - min(modularity, na.rm = TRUE)) / 
                            (max(modularity, na.rm = TRUE) - min(modularity, na.rm = TRUE))),
    stability_norm = ifelse(is.na(stability), 0,
                           (stability - min(stability, na.rm = TRUE)) / 
                           (max(stability, na.rm = TRUE) - min(stability, na.rm = TRUE))),
    singleton_penalty = n_singletons / n_clusters,
    large_cluster_penalty = largest_cluster_prop
  ) %>%
  mutate(
    composite_score = (silhouette_norm + modularity_norm + stability_norm) / 3 - 
                     (singleton_penalty + large_cluster_penalty) / 2
  )

# Find optimal resolution
optimal_idx <- which.max(quality_metrics_norm$composite_score)
optimal_resolution <- quality_metrics_norm$resolution[optimal_idx]

cat("Optimal resolution identified:", optimal_resolution, "\\n")
cat("  - Number of clusters:", quality_metrics$n_clusters[optimal_idx], "\\n")
cat("  - Silhouette score:", round(quality_metrics$silhouette_score[optimal_idx], 3), "\\n")
cat("  - Modularity:", round(quality_metrics$modularity[optimal_idx], 3), "\\n")

# Set optimal clustering as default
optimal_col <- paste0("RNA_snn_res.", optimal_resolution)
if (optimal_col %in% colnames(sobj@meta.data)) {
  Idents(sobj) <- optimal_col
  sobj$optimal_clusters <- sobj@meta.data[[optimal_col]]
}

# ============================================================================
# Cluster Characterization
# ============================================================================

cat("Characterizing optimal clusters...\\n")

# Find markers for optimal clustering
optimal_markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, 
                                 logfc.threshold = 0.25, verbose = FALSE)

# Cluster composition analysis
composition_data <- NULL
if ("celltype" %in% colnames(sobj@meta.data)) {
  composition_data <- analyze_cluster_composition(sobj, "celltype")
} else if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
  composition_data <- analyze_cluster_composition(sobj, "preliminary_celltype")
}

# Create cluster summary
cluster_summary <- sobj@meta.data %>%
  group_by(optimal_clusters) %>%
  summarise(
    n_cells = n(),
    mean_nFeature = mean(nFeature_RNA),
    mean_nCount = mean(nCount_RNA),
    .groups = 'drop'
  )

# Add top markers
if (nrow(optimal_markers) > 0) {
  top_markers_summary <- optimal_markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) %>%
    summarise(top_markers = paste(gene, collapse = ", "), .groups = 'drop')
  
  cluster_summary <- left_join(cluster_summary, top_markers_summary, 
                              by = c("optimal_clusters" = "cluster"))
}

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating comprehensive report...\\n")

# Summary statistics
n_cells <- ncol(sobj)
n_resolutions_tested <- length(RESOLUTION_RANGE)
optimal_n_clusters <- quality_metrics$n_clusters[optimal_idx]

# Generate report
report_text <- paste0(
  "Multi-Resolution Clustering Analysis Report\\n",
  "==========================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Dataset: ", n_cells, " cells\\n\\n",
  "Analysis Parameters:\\n",
  "  - Resolution range: ", min(RESOLUTION_RANGE), " - ", max(RESOLUTION_RANGE), "\\n",
  "  - Resolution step: ", RESOLUTION_RANGE[2] - RESOLUTION_RANGE[1], "\\n",
  "  - PCA dimensions: ", paste(range(PCA_DIMS), collapse = "-"), "\\n",
  "  - Min cluster size threshold: ", MIN_CLUSTER_SIZE, "\\n\\n",
  "Results Summary:\\n",
  "  - Resolutions tested: ", n_resolutions_tested, "\\n",
  "  - Optimal resolution: ", optimal_resolution, "\\n",
  "  - Optimal number of clusters: ", optimal_n_clusters, "\\n",
  "  - Quality metrics calculated: ", paste(METRICS_TO_CALCULATE, collapse = ", "), "\\n\\n",
  "Optimal Resolution Metrics:\\n",
  "  - Silhouette score: ", round(quality_metrics$silhouette_score[optimal_idx], 3), "\\n",
  "  - Modularity: ", round(quality_metrics$modularity[optimal_idx], 3), "\\n",
  "  - Small clusters: ", quality_metrics$n_singletons[optimal_idx], "\\n",
  "  - Largest cluster proportion: ", round(quality_metrics$largest_cluster_prop[optimal_idx], 3), "\\n\\n",
  "Output Files:\\n",
  "  - Quality metrics: multi_resolution_quality_metrics.csv\\n",
  "  - Cluster summary: optimal_clustering_summary.csv\\n",
  "  - Marker genes: optimal_clustering_markers.csv\\n",
  "  - Processed object: multi_resolution_analyzed_object.rds\\n\\n",
  "Visualization Files:\\n",
  "  - multi_resolution_overview.pdf: Resolution comparison\\n",
  "  - clustering_quality_metrics.pdf: Quality metric trends\\n",
  "  - cluster_count_vs_resolution.pdf: Cluster count analysis\\n",
  "  - clustree_resolution_tree.pdf: Resolution hierarchy tree\\n\\n",
  "Recommendations:\\n",
  "1. Use resolution ", optimal_resolution, " for downstream analysis\\n",
  "2. Review small clusters for potential merging\\n",
  "3. Validate cluster assignments with biological knowledge\\n",
  "4. Consider higher resolutions for subclustering specific populations\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "multi_resolution_analysis_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\\n")

# Save quality metrics
write.csv(quality_metrics, file.path(OUTPUT_DIR, "multi_resolution_quality_metrics.csv"), 
          row.names = FALSE)

# Save cluster summary
write.csv(cluster_summary, file.path(OUTPUT_DIR, "optimal_clustering_summary.csv"), 
          row.names = FALSE)

# Save marker genes
if (nrow(optimal_markers) > 0) {
  write.csv(optimal_markers, file.path(OUTPUT_DIR, "optimal_clustering_markers.csv"), 
            row.names = FALSE)
}

# Save composition data
if (!is.null(composition_data)) {
  write.csv(composition_data, file.path(OUTPUT_DIR, "cluster_composition_analysis.csv"), 
            row.names = FALSE)
}

# Save processed Seurat object
saveRDS(sobj, file.path(OUTPUT_DIR, "multi_resolution_analyzed_object.rds"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_multi_resolution_analysis.txt"))

cat("\\nMulti-resolution clustering analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nKey findings:\\n")
cat("- Optimal resolution:", optimal_resolution, "\\n")
cat("- Number of optimal clusters:", optimal_n_clusters, "\\n")
cat("- Quality score:", round(quality_metrics_norm$composite_score[optimal_idx], 3), "\\n")

if (nrow(optimal_markers) > 0) {
  cat("- Marker genes identified:", nrow(optimal_markers), "\\n")
}

cat("\\nNext steps:\\n")
cat("1. Review optimal clustering results\\n")
cat("2. Validate clusters with domain knowledge\\n")
cat("3. Use optimal resolution for downstream analysis\\n")