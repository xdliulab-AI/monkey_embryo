#!/usr/bin/env Rscript

# ============================================================================
# Initial Clustering Analysis for Spatial Transcriptomics Data
# ============================================================================
# 
# This script performs initial clustering analysis on spatial transcriptomics 
# data from primate embryos at Carnegie Stages 9 and 10. It includes quality
# control, normalization, dimensionality reduction, and clustering steps.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
})

# Set parameters
# These can be modified based on your specific dataset
MIN_FEATURES <- 200        # Minimum number of features per spot
MAX_FEATURES <- 8000       # Maximum number of features per spot  
MITOCHONDRIAL_THRESHOLD <- 20  # Maximum percentage of mitochondrial genes
PCA_DIMS <- 1:30          # PCA dimensions to use for clustering
CLUSTERING_RESOLUTION <- 0.5  # Clustering resolution

# Input/Output paths - modify these according to your data structure
INPUT_DATA_PATH <- "../data/spatial_count_matrix.rds"  # Path to count matrix
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories if they don't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Data Loading and Initial Processing
# ============================================================================

cat("Loading spatial transcriptomics data...\n")

# Load the spatial transcriptomics count matrix
# Note: Adjust this section based on your data format
if (file.exists(INPUT_DATA_PATH)) {
  count_matrix <- readRDS(INPUT_DATA_PATH)
} else {
  stop("Input data file not found. Please check the INPUT_DATA_PATH.")
}

# Create Seurat object
sobj <- CreateSeuratObject(
  counts = count_matrix,
  project = "PrimateEmbryo_Spatial",
  min.cells = 3,
  min.features = MIN_FEATURES
)

# Add metadata if available (stage, slice information, etc.)
# Modify this section based on your metadata structure
if (file.exists("../data/metadata.csv")) {
  metadata <- read.csv("../data/metadata.csv", row.names = 1)
  # Ensure metadata rows match Seurat object cells
  common_cells <- intersect(rownames(metadata), colnames(sobj))
  sobj <- sobj[, common_cells]
  sobj <- AddMetaData(sobj, metadata[common_cells, ])
}

cat("Initial data dimensions:", dim(sobj), "\n")

# ============================================================================
# Quality Control
# ============================================================================

cat("Performing quality control analysis...\n")

# Calculate mitochondrial gene percentage
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

# Calculate ribosomal gene percentage (optional)
sobj[["percent.rb"]] <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")

# Visualize QC metrics
qc_plot <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                   ncol = 3, pt.size = 0.1)
ggsave(file.path(PLOTS_DIR, "QC_violin_plot.pdf"), qc_plot, 
       width = 12, height = 6)

# Feature scatter plots
feature_plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
feature_plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc_scatter <- feature_plot1 + feature_plot2
ggsave(file.path(PLOTS_DIR, "QC_scatter_plot.pdf"), qc_scatter, 
       width = 12, height = 6)

# ============================================================================
# Filtering
# ============================================================================

cat("Applying quality control filters...\n")

# Store pre-filtering counts
n_cells_before <- ncol(sobj)
n_features_before <- nrow(sobj)

# Apply filters
sobj <- subset(sobj, subset = nFeature_RNA > MIN_FEATURES & 
                             nFeature_RNA < MAX_FEATURES & 
                             percent.mt < MITOCHONDRIAL_THRESHOLD)

# Report filtering results
n_cells_after <- ncol(sobj)
n_features_after <- nrow(sobj)

cat("Filtering results:\n")
cat("  Cells: ", n_cells_before, " -> ", n_cells_after, 
    " (", round(100 * n_cells_after / n_cells_before, 2), "%)\n")
cat("  Features: ", n_features_before, " -> ", n_features_after, 
    " (", round(100 * n_features_after / n_features_before, 2), "%)\n")

# ============================================================================
# Normalization and Scaling
# ============================================================================

cat("Normalizing and scaling data...\n")

# Normalize the data using SCTransform (recommended for spatial data)
sobj <- SCTransform(sobj, assay = "RNA", verbose = FALSE)

# Alternative: Standard normalization (uncomment if preferred)
# sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
# sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
# sobj <- ScaleData(sobj, features = rownames(sobj))

# ============================================================================
# Dimensionality Reduction
# ============================================================================

cat("Performing dimensionality reduction...\n")

# Principal Component Analysis
sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)

# Visualize PCA results
pca_plot <- DimPlot(sobj, reduction = "pca", group.by = "orig.ident")
ggsave(file.path(PLOTS_DIR, "PCA_plot.pdf"), pca_plot, width = 8, height = 6)

# Elbow plot to determine optimal number of PCs
elbow_plot <- ElbowPlot(sobj, ndims = 50)
ggsave(file.path(PLOTS_DIR, "PCA_elbow_plot.pdf"), elbow_plot, width = 8, height = 6)

# UMAP embedding
sobj <- RunUMAP(sobj, reduction = "pca", dims = PCA_DIMS, verbose = FALSE)

# ============================================================================
# Clustering
# ============================================================================

cat("Performing clustering analysis...\n")

# Find neighbors
sobj <- FindNeighbors(sobj, reduction = "pca", dims = PCA_DIMS, verbose = FALSE)

# Find clusters at specified resolution
sobj <- FindClusters(sobj, resolution = CLUSTERING_RESOLUTION, verbose = FALSE)

# Visualize clustering results
cluster_umap <- DimPlot(sobj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle(paste("UMAP Clustering (resolution =", CLUSTERING_RESOLUTION, ")"))
ggsave(file.path(PLOTS_DIR, paste0("UMAP_clustering_res", CLUSTERING_RESOLUTION, ".pdf")), 
       cluster_umap, width = 10, height = 8)

# If spatial coordinates are available, create spatial plot
if ("Spatial" %in% names(sobj@reductions) || 
    any(c("x", "y", "imagerow", "imagecol") %in% colnames(sobj@meta.data))) {
  
  # Create spatial plot (adjust based on your coordinate column names)
  if ("Spatial" %in% names(sobj@reductions)) {
    spatial_plot <- DimPlot(sobj, reduction = "Spatial", label = TRUE, pt.size = 0.3) +
      ggtitle("Spatial Clustering")
  } else {
    # If coordinates are in metadata, create custom spatial plot
    spatial_data <- data.frame(
      x = sobj@meta.data$x,  # Adjust column name as needed
      y = sobj@meta.data$y,  # Adjust column name as needed
      cluster = Idents(sobj)
    )
    spatial_plot <- ggplot(spatial_data, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.3) +
      coord_fixed() +
      theme_minimal() +
      ggtitle("Spatial Clustering")
  }
  
  ggsave(file.path(PLOTS_DIR, "Spatial_clustering.pdf"), 
         spatial_plot, width = 12, height = 10)
}

# ============================================================================
# Cluster Statistics and Summary
# ============================================================================

cat("Generating cluster statistics...\n")

# Number of clusters
n_clusters <- length(unique(Idents(sobj)))
cat("Number of clusters identified:", n_clusters, "\n")

# Cluster sizes
cluster_sizes <- table(Idents(sobj))
cat("Cluster sizes:\n")
print(cluster_sizes)

# Create summary table
cluster_summary <- data.frame(
  Cluster = names(cluster_sizes),
  n_cells = as.numeric(cluster_sizes),
  percentage = round(100 * as.numeric(cluster_sizes) / ncol(sobj), 2)
)

write.csv(cluster_summary, file.path(OUTPUT_DIR, "cluster_summary.csv"), 
          row.names = FALSE)

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\n")

# Save the processed Seurat object
saveRDS(sobj, file.path(OUTPUT_DIR, "sobj_initial_clustering.rds"))

# Save cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(sobj),
  cluster = as.character(Idents(sobj)),
  stringsAsFactors = FALSE
)
write.csv(cluster_assignments, file.path(OUTPUT_DIR, "cluster_assignments.csv"), 
          row.names = FALSE)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_initial_clustering.txt"))

cat("Initial clustering analysis completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

# ============================================================================
# Generate Summary Report
# ============================================================================

summary_text <- paste0(
  "Initial Clustering Analysis Summary\n",
  "===================================\n\n",
  "Input data: ", INPUT_DATA_PATH, "\n",
  "Date: ", Sys.Date(), "\n\n",
  "Data dimensions after filtering:\n",
  "  - Cells: ", ncol(sobj), "\n",
  "  - Features: ", nrow(sobj), "\n\n",
  "Clustering parameters:\n",
  "  - PCA dimensions: ", paste(range(PCA_DIMS), collapse = "-"), "\n",
  "  - Resolution: ", CLUSTERING_RESOLUTION, "\n",
  "  - Number of clusters: ", n_clusters, "\n\n",
  "Output files:\n",
  "  - Seurat object: sobj_initial_clustering.rds\n",
  "  - Cluster assignments: cluster_assignments.csv\n",
  "  - Cluster summary: cluster_summary.csv\n",
  "  - Plots: QC and clustering visualization PDFs\n"
)

writeLines(summary_text, file.path(OUTPUT_DIR, "analysis_summary.txt"))

cat("\nAnalysis summary saved to: analysis_summary.txt\n")