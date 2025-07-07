#!/usr/bin/env Rscript

# ============================================================================
# Harmony Integration Analysis for Multi-Sample Spatial Data
# ============================================================================
# 
# This script performs batch correction and integration of multiple spatial
# transcriptomics samples using Harmony. It handles spatial coordinate
# transformations and creates unified embeddings for comparative analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
  library(harmony)
  library(tidyverse)
  library(dplyr)
  library(Matrix)
  library(patchwork)
  library(loupeR)
})

# Parameters
INTEGRATION_METHOD <- "HarmonyIntegration"  # Integration method
IDENTIFIER <- "harmony"                     # Method identifier
PCA_DIMS <- 1:30                           # PCA dimensions
CLUSTERING_RESOLUTION <- seq(0.1, 2, 0.1)  # Multiple clustering resolutions
GROUP_GAP <- 100                           # Gap between developmental stages

# Sample information
CS9_SAMPLES <- c("D03251D414", "D03251E211", "D03251G212", "D03251G513", 
                 "D03254F211", "D03257F312", "D03257G311", "D03257G413", 
                 "D03259G313", "D03259G613")

CS10_SAMPLES <- c("D03251F413", "D03251G111", "D03252G213", "D03254C314", 
                  "D03255E612", "D03257C314", "D03257C412", "D03257D411", 
                  "D03257E112", "D03257G212")

ALL_SAMPLES <- c(CS9_SAMPLES, CS10_SAMPLES)

# Input/Output paths
INPUT_DIR <- "../data"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Data Loading and Spatial Transformation
# ============================================================================

cat("Loading and transforming spatial data...\\n")

# Function to load and transform individual samples
load_and_transform_sample <- function(sample_id) {
  file_path <- file.path(INPUT_DIR, paste0("sobj_", sample_id, ".rds"))
  
  if (!file.exists(file_path)) {
    cat("Warning: File not found for sample", sample_id, "\\n")
    return(NULL)
  }
  
  tryCatch({
    # Load sample
    sobj <- readRDS(file_path)
    
    # Add Carnegie Stage information
    sobj$CS <- ifelse(sample_id %in% CS9_SAMPLES, "CS9", "CS10")
    sobj$sample_id <- sample_id
    
    # Extract spatial coordinates
    spatial_coords <- sobj@reductions[["Spatial"]]@cell.embeddings
    spatial_coords <- as.data.frame(spatial_coords)
    spatial_coords$sample_id <- sample_id
    
    # Apply transformation for CS9 samples (rotation and flip)
    if (sample_id %in% CS9_SAMPLES) {
      spatial_coords <- spatial_coords %>%
        mutate(
          # Rotation: new X = -Y, new Y = X
          rotated_Spatial_1 = -Spatial_2,
          rotated_Spatial_2 = Spatial_1
        ) %>%
        mutate(
          # Flip Y-axis
          flipped_Spatial_2 = -rotated_Spatial_2
        ) %>%
        mutate(
          # Update coordinates
          Spatial_1 = rotated_Spatial_1,
          Spatial_2 = flipped_Spatial_2
        ) %>%
        select(Spatial_1, Spatial_2, sample_id)
    }
    
    return(list(sobj = sobj, coords = spatial_coords))
    
  }, error = function(e) {
    cat("Error loading sample", sample_id, ":", e$message, "\\n")
    return(NULL)
  })
}

# Load all samples
cat("Loading", length(ALL_SAMPLES), "samples...\\n")
sample_data <- map(ALL_SAMPLES, load_and_transform_sample)
names(sample_data) <- ALL_SAMPLES

# Remove NULL entries (failed loads)
sample_data <- sample_data[!sapply(sample_data, is.null)]
loaded_samples <- names(sample_data)

cat("Successfully loaded", length(loaded_samples), "samples\\n")

# Extract Seurat objects and spatial coordinates
sobj_list <- map(sample_data, ~ .x$sobj)
spatial_coords_list <- map(sample_data, ~ .x$coords)

# Combine spatial coordinates
combined_spatial <- bind_rows(spatial_coords_list)

# ============================================================================
# Spatial Coordinate Transformations
# ============================================================================

cat("Creating spatial coordinate transformations...\\n")

# 1. SpatialCenter - Center each slice at (0,0)
spatial_center <- combined_spatial %>%
  group_by(sample_id) %>%
  mutate(
    mean_Spatial_1 = mean(Spatial_1),
    mean_Spatial_2 = mean(Spatial_2),
    adjusted_Spatial_1 = Spatial_1 - mean_Spatial_1,
    adjusted_Spatial_2 = Spatial_2 - mean_Spatial_2
  ) %>%
  ungroup() %>%
  select(-mean_Spatial_1, -mean_Spatial_2)

# 2. SpatialCenterByCS - Center CS9 and CS10 groups separately
cs9_coords <- spatial_center %>% filter(sample_id %in% CS9_SAMPLES)
cs10_coords <- spatial_center %>% filter(sample_id %in% CS10_SAMPLES)

# Center-align by group
cs9_coords <- cs9_coords %>%
  mutate(
    group_mean_Spatial_1 = mean(adjusted_Spatial_1),
    group_mean_Spatial_2 = mean(adjusted_Spatial_2),
    group_adjusted_Spatial_1 = adjusted_Spatial_1 - group_mean_Spatial_1 - GROUP_GAP/2,
    group_adjusted_Spatial_2 = adjusted_Spatial_2 - group_mean_Spatial_2
  ) %>%
  select(-group_mean_Spatial_1, -group_mean_Spatial_2)

cs10_coords <- cs10_coords %>%
  mutate(
    group_mean_Spatial_1 = mean(adjusted_Spatial_1),
    group_mean_Spatial_2 = mean(adjusted_Spatial_2),
    group_adjusted_Spatial_1 = adjusted_Spatial_1 - group_mean_Spatial_1 + GROUP_GAP/2,
    group_adjusted_Spatial_2 = adjusted_Spatial_2 - group_mean_Spatial_2
  ) %>%
  select(-group_mean_Spatial_1, -group_mean_Spatial_2)

spatial_center_by_cs <- bind_rows(cs9_coords, cs10_coords)

# 3. Spatial45 - 4x5 grid layout
column_count <- 5
row_count <- 4
custom_vertical_shifts <- c(0.9, 1.2, 1.2)

horizontal_shift <- (max(spatial_center$adjusted_Spatial_1) - 
                    min(spatial_center$adjusted_Spatial_1)) * 1.5
base_vertical_shift <- max(spatial_center$adjusted_Spatial_2) - 
                      min(spatial_center$adjusted_Spatial_2)

spatial_45 <- spatial_center
unique_ids <- unique(spatial_center$sample_id)

for (i in seq_along(unique_ids)) {
  row <- (i - 1) %/% column_count
  col <- (i - 1) %% column_count
  
  vertical_shift <- ifelse(row == 0, 0, 
                          sum(custom_vertical_shifts[1:row]) * base_vertical_shift)
  
  sample_mask <- spatial_45$sample_id == unique_ids[i]
  spatial_45$Spatial_1[sample_mask] <- spatial_45$adjusted_Spatial_1[sample_mask] + 
                                      col * horizontal_shift
  spatial_45$Spatial_2[sample_mask] <- spatial_45$adjusted_Spatial_2[sample_mask] - 
                                      vertical_shift
}

# ============================================================================
# Merge Seurat Objects
# ============================================================================

cat("Merging Seurat objects...\\n")

# Merge all objects
merged_sobj <- Reduce(function(x, y) merge(x, y), sobj_list)

# Join layers and clean up
merged_sobj[["RNA"]] <- JoinLayers(merged_sobj[["RNA"]])
merged_sobj@assays[["RNA"]]@layers[["scale.data"]] <- NULL
merged_sobj@assays[["RNA"]]@layers[["data"]] <- NULL

# Clean metadata
md <- merged_sobj@meta.data
md <- md[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "CS", "sample_id")]
merged_sobj@meta.data <- md

# Split assay by sample
merged_sobj[["RNA"]] <- split(merged_sobj[["RNA"]], f = merged_sobj$sample_id)

# ============================================================================
# Preprocessing and Integration
# ============================================================================

cat("Preprocessing for integration...\\n")

# Standard preprocessing
merged_sobj <- NormalizeData(merged_sobj, verbose = FALSE)
merged_sobj <- FindVariableFeatures(merged_sobj, verbose = FALSE)
merged_sobj <- ScaleData(merged_sobj, verbose = FALSE)
merged_sobj <- RunPCA(merged_sobj, verbose = FALSE)

# Harmony integration
cat("Performing Harmony integration...\\n")

reduction_name <- paste('integrated', IDENTIFIER, sep = '.')

merged_sobj <- IntegrateLayers(
  object = merged_sobj, 
  method = INTEGRATION_METHOD,
  new.reduction = reduction_name,
  verbose = FALSE
)

# Re-join layers
merged_sobj[["RNA"]] <- JoinLayers(merged_sobj[["RNA"]])

# ============================================================================
# Downstream Analysis
# ============================================================================

cat("Performing downstream analysis...\\n")

# Clustering and UMAP
merged_sobj <- FindNeighbors(merged_sobj, reduction = reduction_name, 
                            dims = PCA_DIMS, verbose = FALSE)
merged_sobj <- FindClusters(merged_sobj, resolution = CLUSTERING_RESOLUTION, 
                           verbose = FALSE)
merged_sobj <- RunUMAP(merged_sobj, dims = PCA_DIMS, reduction = reduction_name, 
                      verbose = FALSE)

# ============================================================================
# Add Spatial Reductions
# ============================================================================

cat("Adding spatial coordinate reductions...\\n")

# Function to add spatial reduction
add_spatial_reduction <- function(sobj, coords_df, reduction_name) {
  # Prepare coordinate matrix
  coord_matrix <- as.matrix(coords_df[, c("Spatial_1", "Spatial_2")])
  rownames(coord_matrix) <- colnames(sobj)
  
  # Create reduction object
  spatial_reduction <- CreateDimReducObject(
    embeddings = coord_matrix,
    key = "Spatial_",
    assay = DefaultAssay(sobj)
  )
  
  # Add to Seurat object
  sobj@reductions[[reduction_name]] <- spatial_reduction
  
  return(sobj)
}

# Add all spatial coordinate systems
merged_sobj <- add_spatial_reduction(merged_sobj, 
                                    spatial_center[, c("Spatial_1", "Spatial_2")] %>%
                                      rename(Spatial_1 = adjusted_Spatial_1, 
                                             Spatial_2 = adjusted_Spatial_2),
                                    "SpatialCenter")

merged_sobj <- add_spatial_reduction(merged_sobj, 
                                    spatial_center_by_cs[, c("Spatial_1", "Spatial_2")] %>%
                                      rename(Spatial_1 = group_adjusted_Spatial_1, 
                                             Spatial_2 = group_adjusted_Spatial_2),
                                    "SpatialCenterByCS")

merged_sobj <- add_spatial_reduction(merged_sobj, 
                                    spatial_45[, c("Spatial_1", "Spatial_2")],
                                    "Spatial45")

# ============================================================================
# Visualization
# ============================================================================

cat("Creating visualizations...\\n")

# UMAP plots
umap_by_sample <- DimPlot(merged_sobj, reduction = "umap", 
                         group.by = "sample_id", pt.size = 0.3) +
  ggtitle("Harmony Integration by Sample")

umap_by_cs <- DimPlot(merged_sobj, reduction = "umap", 
                     group.by = "CS", pt.size = 0.3) +
  ggtitle("Harmony Integration by Carnegie Stage")

umap_by_cluster <- DimPlot(merged_sobj, reduction = "umap", 
                          group.by = "RNA_snn_res.0.5", 
                          label = TRUE, pt.size = 0.3) +
  ggtitle("Harmony Integration Clusters")

# Spatial plots
spatial_by_sample <- DimPlot(merged_sobj, reduction = "SpatialCenter", 
                            group.by = "sample_id", pt.size = 0.3) +
  ggtitle("Spatial Distribution by Sample") +
  coord_fixed()

spatial_by_cs <- DimPlot(merged_sobj, reduction = "SpatialCenterByCS", 
                        group.by = "CS", pt.size = 0.3) +
  ggtitle("Spatial Distribution by Carnegie Stage") +
  coord_fixed()

spatial_grid <- DimPlot(merged_sobj, reduction = "Spatial45", 
                       group.by = "sample_id", pt.size = 0.3) +
  ggtitle("4x5 Grid Layout") +
  coord_fixed()

# Save plots
ggsave(file.path(PLOTS_DIR, "harmony_integration_umap.pdf"), 
       (umap_by_sample | umap_by_cs) / umap_by_cluster, 
       width = 16, height = 12)

ggsave(file.path(PLOTS_DIR, "harmony_integration_spatial.pdf"), 
       (spatial_by_sample | spatial_by_cs) / spatial_grid, 
       width = 20, height = 16)

# ============================================================================
# Quality Assessment
# ============================================================================

cat("Assessing integration quality...\\n")

# Calculate integration metrics
integration_metrics <- data.frame(
  metric = c("Total cells", "Total samples", "CS9 samples", "CS10 samples", 
             "Clusters (res 0.5)", "Features"),
  value = c(ncol(merged_sobj), length(loaded_samples), 
            length(CS9_SAMPLES), length(CS10_SAMPLES),
            length(unique(merged_sobj$RNA_snn_res.0.5)), 
            nrow(merged_sobj))
)

# Sample composition per cluster
cluster_composition <- merged_sobj@meta.data %>%
  group_by(RNA_snn_res.0.5, CS) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = CS, values_from = n_cells, values_fill = 0) %>%
  mutate(
    total_cells = CS9 + CS10,
    cs9_prop = CS9 / total_cells,
    cs10_prop = CS10 / total_cells,
    balance_score = 1 - abs(cs9_prop - cs10_prop)
  )

# Save quality metrics
write.csv(integration_metrics, file.path(OUTPUT_DIR, "harmony_integration_metrics.csv"), 
          row.names = FALSE)
write.csv(cluster_composition, file.path(OUTPUT_DIR, "harmony_cluster_composition.csv"), 
          row.names = FALSE)

# ============================================================================
# Generate Loupe File
# ============================================================================

cat("Generating Loupe file...\\n")

# Generate unique barcodes
n_cells <- ncol(merged_sobj)
barcodes_file <- file.path(INPUT_DIR, "barcodes.tsv")

if (file.exists(barcodes_file)) {
  all_barcodes <- readr::read_tsv(barcodes_file, col_names = FALSE)[[1]]
  set.seed(123)
  unique_barcodes <- sample(all_barcodes, size = n_cells, replace = FALSE)
  colnames(merged_sobj) <- unique_barcodes
}

# Set factor levels
merged_sobj$sample_id <- factor(merged_sobj$sample_id, levels = ALL_SAMPLES)
merged_sobj$CS <- factor(merged_sobj$CS, levels = c("CS9", "CS10"))

# Create Loupe file
loupe_output <- file.path(OUTPUT_DIR, paste0("loupe_harmony_integration_", 
                                            Sys.Date()))
create_loupe_from_seurat(merged_sobj, output_name = loupe_output, force = TRUE)

# ============================================================================
# Generate Integration Report
# ============================================================================

cat("Generating integration report...\\n")

# Calculate summary statistics
n_total_cells <- ncol(merged_sobj)
n_cs9_cells <- sum(merged_sobj$CS == "CS9")
n_cs10_cells <- sum(merged_sobj$CS == "CS10")
n_clusters <- length(unique(merged_sobj$RNA_snn_res.0.5))
avg_balance_score <- mean(cluster_composition$balance_score, na.rm = TRUE)

# Create report
report_text <- paste0(
  "Harmony Integration Analysis Report\\n",
  "=================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n\\n",
  "Data Summary:\\n",
  "  - Total cells: ", n_total_cells, "\\n",
  "  - CS9 cells: ", n_cs9_cells, " (", round(100 * n_cs9_cells / n_total_cells, 1), "%)\\n",
  "  - CS10 cells: ", n_cs10_cells, " (", round(100 * n_cs10_cells / n_total_cells, 1), "%)\\n",
  "  - Samples loaded: ", length(loaded_samples), "/", length(ALL_SAMPLES), "\\n",
  "  - Clusters identified: ", n_clusters, "\\n\\n",
  "Integration Parameters:\\n",
  "  - Method: ", INTEGRATION_METHOD, "\\n",
  "  - PCA dimensions: ", paste(range(PCA_DIMS), collapse = "-"), "\\n",
  "  - Clustering resolutions: ", paste(range(CLUSTERING_RESOLUTION), collapse = "-"), "\\n\\n",
  "Spatial Transformations:\\n",
  "  - SpatialCenter: Individual slice centering\\n",
  "  - SpatialCenterByCS: Group-wise centering with gap\\n",
  "  - Spatial45: 4x5 grid layout\\n\\n",
  "Integration Quality:\\n",
  "  - Average balance score: ", round(avg_balance_score, 3), "\\n",
  "  - Well-balanced clusters (>0.7): ", sum(cluster_composition$balance_score > 0.7, na.rm = TRUE), 
    "/", nrow(cluster_composition), "\\n\\n",
  "Output Files:\\n",
  "  - Integrated object: harmony_integrated_object.rds\\n",
  "  - Integration metrics: harmony_integration_metrics.csv\\n",
  "  - Cluster composition: harmony_cluster_composition.csv\\n",
  "  - Loupe file: ", basename(loupe_output), "\\n\\n",
  "Visualization Files:\\n",
  "  - harmony_integration_umap.pdf: UMAP visualizations\\n",
  "  - harmony_integration_spatial.pdf: Spatial visualizations\\n\\n",
  "Recommendations:\\n",
  "1. Review cluster composition for balanced representation\\n",
  "2. Validate spatial transformations are appropriate\\n",
  "3. Consider additional batch correction if needed\\n",
  "4. Use integrated data for downstream analysis\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "harmony_integration_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\\n")

# Save integrated object
saveRDS(merged_sobj, file.path(OUTPUT_DIR, "harmony_integrated_object.rds"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_harmony_integration.txt"))

cat("\\nHarmony integration analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nIntegration summary:\\n")
print(integration_metrics)
cat("\\nNext steps:\\n")
cat("1. Review integration quality metrics\\n")
cat("2. Validate spatial coordinate transformations\\n")
cat("3. Proceed with downstream analysis using integrated data\\n")