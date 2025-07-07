#!/usr/bin/env Rscript

# ============================================================================
# scRNA-seq and Spatial Transcriptomics Integration Analysis
# ============================================================================
# 
# This script integrates single-cell RNA-seq data with spatial transcriptomics
# data from primate embryos. It performs data merging, batch correction, and
# creates unified embeddings for downstream analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
  library(tidyverse)
  library(dplyr)
  library(Matrix)
  library(patchwork)
})

# Parameters
MIN_CELLS <- 3            # Minimum cells per feature
MIN_FEATURES <- 200       # Minimum features per cell
PCA_DIMS <- 1:30         # PCA dimensions for integration
CLUSTERING_RESOLUTION <- 0.5  # Clustering resolution
INTEGRATION_METHOD <- "RPCAIntegration"  # Integration method

# Input/Output paths
INPUT_SCRNA_PATH <- "../data/scrna_seurat_object.rds"
INPUT_SPATIAL_PATH <- "../data/spatial_seurat_object.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load and Prepare Data
# ============================================================================

cat("Loading scRNA-seq and spatial transcriptomics data...\\n")

# Load scRNA-seq data
if (file.exists(INPUT_SCRNA_PATH)) {
  scrna_obj <- readRDS(INPUT_SCRNA_PATH)
  cat("Loaded scRNA-seq data with", ncol(scrna_obj), "cells\\n")
} else {
  stop("scRNA-seq data file not found: ", INPUT_SCRNA_PATH)
}

# Load spatial data
if (file.exists(INPUT_SPATIAL_PATH)) {
  spatial_obj <- readRDS(INPUT_SPATIAL_PATH)
  cat("Loaded spatial data with", ncol(spatial_obj), "spots\\n")
} else {
  stop("Spatial data file not found: ", INPUT_SPATIAL_PATH)
}

# Add data type labels
scrna_obj$data_type <- "scRNA-seq"
spatial_obj$data_type <- "Spatial"

# ============================================================================
# Data Preprocessing
# ============================================================================

cat("Preprocessing data for integration...\\n")

# Function to preprocess Seurat objects
preprocess_seurat <- function(sobj, sample_name) {
  # Add sample identity
  sobj$sample_id <- sample_name
  
  # Quality control metrics
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.rb"]] <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
  
  # Filter low quality cells/spots
  sobj <- subset(sobj, subset = nFeature_RNA > MIN_FEATURES & percent.mt < 20)
  
  # Normalize and find variable features
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", 
                              nfeatures = 2000, verbose = FALSE)
  
  return(sobj)
}

# Preprocess both datasets
scrna_obj <- preprocess_seurat(scrna_obj, "scRNA")
spatial_obj <- preprocess_seurat(spatial_obj, "Spatial")

# ============================================================================
# Integration Analysis
# ============================================================================

cat("Performing integration analysis...\\n")

# Create list of objects to integrate
integration_list <- list(scrna = scrna_obj, spatial = spatial_obj)

# Find integration anchors
integration_anchors <- FindIntegrationAnchors(
  object.list = integration_list,
  dims = PCA_DIMS,
  verbose = FALSE
)

# Integrate data
integrated_obj <- IntegrateData(
  anchorset = integration_anchors,
  dims = PCA_DIMS,
  verbose = FALSE
)

# Switch to integrated assay
DefaultAssay(integrated_obj) <- "integrated"

# ============================================================================
# Dimensionality Reduction and Clustering
# ============================================================================

cat("Performing dimensionality reduction and clustering...\\n")

# Scale data
integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)

# PCA
integrated_obj <- RunPCA(integrated_obj, verbose = FALSE)

# UMAP
integrated_obj <- RunUMAP(integrated_obj, dims = PCA_DIMS, verbose = FALSE)

# Clustering
integrated_obj <- FindNeighbors(integrated_obj, dims = PCA_DIMS, verbose = FALSE)
integrated_obj <- FindClusters(integrated_obj, resolution = CLUSTERING_RESOLUTION, 
                              verbose = FALSE)

# ============================================================================
# Visualization
# ============================================================================

cat("Creating integration visualizations...\\n")

# UMAP colored by data type
data_type_plot <- DimPlot(integrated_obj, group.by = "data_type", 
                         pt.size = 0.3, cols = c("red", "blue")) +
  ggtitle("Integration by Data Type")

# UMAP colored by clusters
cluster_plot <- DimPlot(integrated_obj, label = TRUE, pt.size = 0.3) +
  ggtitle("Integrated Clusters")

# UMAP colored by sample
sample_plot <- DimPlot(integrated_obj, group.by = "sample_id", pt.size = 0.3) +
  ggtitle("Integration by Sample")

# Combine plots
integration_overview <- (data_type_plot | cluster_plot) / sample_plot
ggsave(file.path(PLOTS_DIR, "integration_overview.pdf"), 
       integration_overview, width = 16, height = 12)

# Split by data type
split_plot <- DimPlot(integrated_obj, split.by = "data_type", pt.size = 0.3) +
  ggtitle("Integration Split by Data Type")
ggsave(file.path(PLOTS_DIR, "integration_split_by_type.pdf"), 
       split_plot, width = 16, height = 8)

# ============================================================================
# Quality Assessment
# ============================================================================

cat("Assessing integration quality...\\n")

# Calculate integration metrics
integration_metrics <- data.frame(
  data_type = c("scRNA-seq", "Spatial"),
  n_cells = c(
    sum(integrated_obj$data_type == "scRNA-seq"),
    sum(integrated_obj$data_type == "Spatial")
  ),
  n_features = c(
    sum(rowSums(GetAssayData(scrna_obj, slot = "counts")) > 0),
    sum(rowSums(GetAssayData(spatial_obj, slot = "counts")) > 0)
  )
)

# Calculate mixing metrics per cluster
cluster_composition <- integrated_obj@meta.data %>%
  group_by(seurat_clusters, data_type) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = data_type, values_from = n_cells, values_fill = 0) %>%
  mutate(
    total_cells = `scRNA-seq` + Spatial,
    scrna_prop = `scRNA-seq` / total_cells,
    spatial_prop = Spatial / total_cells,
    mixing_score = 1 - abs(scrna_prop - spatial_prop)
  )

# Save integration metrics
write.csv(integration_metrics, file.path(OUTPUT_DIR, "integration_metrics.csv"), 
          row.names = FALSE)
write.csv(cluster_composition, file.path(OUTPUT_DIR, "cluster_composition.csv"), 
          row.names = FALSE)

# ============================================================================
# Marker Gene Analysis
# ============================================================================

cat("Identifying integration-specific markers...\\n")

# Find markers that distinguish data types
Idents(integrated_obj) <- "data_type"
data_type_markers <- FindAllMarkers(integrated_obj, only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.25, 
                                   verbose = FALSE)

# Find cluster markers
Idents(integrated_obj) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(integrated_obj, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.25, 
                                 verbose = FALSE)

# Save marker results
write.csv(data_type_markers, file.path(OUTPUT_DIR, "data_type_markers.csv"), 
          row.names = FALSE)
write.csv(cluster_markers, file.path(OUTPUT_DIR, "integrated_cluster_markers.csv"), 
          row.names = FALSE)

# ============================================================================
# Cell Type Transfer
# ============================================================================

cat("Performing cell type label transfer...\\n")

# If cell type annotations exist in scRNA-seq data
if ("celltype" %in% colnames(scrna_obj@meta.data)) {
  # Find transfer anchors
  transfer_anchors <- FindTransferAnchors(
    reference = scrna_obj,
    query = spatial_obj,
    dims = PCA_DIMS,
    verbose = FALSE
  )
  
  # Transfer cell type labels
  predictions <- TransferData(
    anchorset = transfer_anchors,
    refdata = scrna_obj$celltype,
    dims = PCA_DIMS,
    verbose = FALSE
  )
  
  # Add predictions to spatial object
  spatial_obj <- AddMetaData(spatial_obj, predictions)
  
  # Add transferred labels to integrated object
  transfer_results <- data.frame(
    cell_id = colnames(spatial_obj),
    predicted_celltype = spatial_obj$predicted.id,
    prediction_score = spatial_obj$prediction.score.max
  )
  
  # Map to integrated object
  integrated_obj$predicted_celltype <- NA
  spatial_cells <- integrated_obj$data_type == "Spatial"
  integrated_obj$predicted_celltype[spatial_cells] <- 
    transfer_results$predicted_celltype[match(
      colnames(integrated_obj)[spatial_cells], 
      transfer_results$cell_id
    )]
  
  # Copy original annotations for scRNA-seq cells
  scrna_cells <- integrated_obj$data_type == "scRNA-seq"
  if ("celltype" %in% colnames(integrated_obj@meta.data)) {
    integrated_obj$predicted_celltype[scrna_cells] <- 
      integrated_obj$celltype[scrna_cells]
  }
  
  # Create transfer visualization
  transfer_plot <- DimPlot(integrated_obj, group.by = "predicted_celltype", 
                          label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle("Transferred Cell Type Labels")
  
  ggsave(file.path(PLOTS_DIR, "cell_type_transfer.pdf"), 
         transfer_plot, width = 14, height = 10)
  
  # Save transfer results
  write.csv(transfer_results, file.path(OUTPUT_DIR, "cell_type_transfer_results.csv"), 
            row.names = FALSE)
}

# ============================================================================
# Generate Integration Report
# ============================================================================

cat("Generating integration report...\\n")

# Calculate summary statistics
n_total_cells <- ncol(integrated_obj)
n_scrna_cells <- sum(integrated_obj$data_type == "scRNA-seq")
n_spatial_cells <- sum(integrated_obj$data_type == "Spatial")
n_clusters <- length(unique(integrated_obj$seurat_clusters))

# Create report
report_text <- paste0(
  "scRNA-seq and Spatial Integration Analysis Report\\n",
  "==============================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n\\n",
  "Data Summary:\\n",
  "  - Total cells/spots: ", n_total_cells, "\\n",
  "  - scRNA-seq cells: ", n_scrna_cells, " (", 
    round(100 * n_scrna_cells / n_total_cells, 1), "%)\\n",
  "  - Spatial spots: ", n_spatial_cells, " (", 
    round(100 * n_spatial_cells / n_total_cells, 1), "%)\\n",
  "  - Integrated clusters: ", n_clusters, "\\n\\n",
  "Integration Parameters:\\n",
  "  - Integration method: ", INTEGRATION_METHOD, "\\n",
  "  - PCA dimensions: ", paste(range(PCA_DIMS), collapse = "-"), "\\n",
  "  - Clustering resolution: ", CLUSTERING_RESOLUTION, "\\n\\n",
  "Output Files:\\n",
  "  - Integrated object: integrated_scrna_spatial.rds\\n",
  "  - Integration metrics: integration_metrics.csv\\n",
  "  - Cluster composition: cluster_composition.csv\\n",
  "  - Data type markers: data_type_markers.csv\\n",
  "  - Cluster markers: integrated_cluster_markers.csv\\n"
)

if ("celltype" %in% colnames(scrna_obj@meta.data)) {
  report_text <- paste0(report_text,
    "  - Cell type transfer: cell_type_transfer_results.csv\\n"
  )
}

report_text <- paste0(report_text,
  "\\nVisualization Files:\\n",
  "  - integration_overview.pdf: Main integration visualization\\n",
  "  - integration_split_by_type.pdf: Split view by data type\\n"
)

if ("celltype" %in% colnames(scrna_obj@meta.data)) {
  report_text <- paste0(report_text,
    "  - cell_type_transfer.pdf: Transferred cell type labels\\n"
  )
}

report_text <- paste0(report_text,
  "\\nIntegration Quality:\\n",
  "  - Average mixing score: ", round(mean(cluster_composition$mixing_score, na.rm = TRUE), 3), "\\n",
  "  - Well-mixed clusters (>0.7): ", sum(cluster_composition$mixing_score > 0.7, na.rm = TRUE), "/", nrow(cluster_composition), "\\n\\n",
  "Recommendations:\\n",
  "1. Review cluster composition for balanced integration\\n",
  "2. Validate transferred cell types with known markers\\n",
  "3. Consider batch effects if integration quality is poor\\n",
  "4. Use integrated data for downstream comparative analysis\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "integration_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving integrated results...\\n")

# Save integrated object
saveRDS(integrated_obj, file.path(OUTPUT_DIR, "integrated_scrna_spatial.rds"))

# Save individual processed objects
saveRDS(scrna_obj, file.path(OUTPUT_DIR, "processed_scrna.rds"))
saveRDS(spatial_obj, file.path(OUTPUT_DIR, "processed_spatial.rds"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_integration.txt"))

cat("\\nIntegration analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nIntegration quality metrics:\\n")
print(integration_metrics)
cat("\\nNext steps:\\n")
cat("1. Review integration quality and cluster composition\\n")
cat("2. Validate cell type transfers using known markers\\n")
cat("3. Use integrated data for comparative analysis\\n")