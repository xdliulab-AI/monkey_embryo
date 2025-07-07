#!/usr/bin/env Rscript

# ============================================================================
# Lineage-Specific Subclustering Analysis
# ============================================================================
# 
# This script performs subclustering analysis on specific cell lineages or
# cell types of interest. It extracts cells from integrated datasets and
# performs refined clustering to identify subtypes and states.
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
  library(loupeR)
  library(patchwork)
})

# Parameters
INTEGRATION_METHOD <- "RPCAIntegration"  # Integration method: RPCAIntegration, HarmonyIntegration, FastMNNIntegration
IDENTIFIER <- "rpca"                     # Method identifier for file naming
PCA_DIMS <- 1:30                        # PCA dimensions to use
CLUSTERING_RESOLUTIONS <- seq(0.1, 3.0, 0.1)  # Clustering resolutions to test
MIN_CELLS_PER_SAMPLE <- 20              # Minimum cells per sample to include

# Cell type selection options - modify based on your analysis needs
LINEAGE_DEFINITIONS <- list(
  "Neural_lineage" = c("Brain_and_spinal_cord", "Neural_ectoderm", "Floor_plate", 
                       "Neural_crest", "CNS_neurons", "PNS_neurons"),
  "Mesodermal_lineage" = c("Paraxial_mesoderm", "Caudal_mesoderm", "Lateral_plate_mesoderm",
                          "Heart", "Cardiac_progenitors", "Endothelial_cells", "Blood_progenitor"),
  "Endodermal_lineage" = c("Gut_tube", "Foregut", "Midgut", "Hindgut", "Yolk_sac_endoderm"),
  "Axial_structures" = c("Notochord", "Prechordal_plate"),
  "Extraembryonic" = c("Allantois", "Amnion", "Yolk_sac_mesoderm", "Yolk_sac_endoderm")
)

# Input/Output paths
INPUT_SEURAT <- "../data/integrated_annotated_object.rds"
LINEAGE_TYPE <- "Neural_lineage"  # Change this to select different lineages
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Parse command line arguments if provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  LINEAGE_TYPE <- args[1]
}

if (length(args) > 1) {
  INTEGRATION_METHOD <- args[2]
  IDENTIFIER <- tolower(gsub("Integration", "", INTEGRATION_METHOD))
}

# ============================================================================
# Load Data and Cell Selection
# ============================================================================

cat("Loading integrated Seurat object...\\n")

if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Determine cell type column
cell_type_col <- "celltype"
if (!"celltype" %in% colnames(sobj@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "preliminary_celltype"
  } else if ("cell_type" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "cell_type"
  } else {
    stop("No suitable cell type annotation found")
  }
}

cat("Using cell type annotation column:", cell_type_col, "\\n")

# Select cell types for the specified lineage
if (!LINEAGE_TYPE %in% names(LINEAGE_DEFINITIONS)) {
  # If lineage type is not predefined, treat as comma-separated cell types
  selected_celltypes <- unlist(strsplit(LINEAGE_TYPE, ","))
  cat("Using custom cell type selection:", paste(selected_celltypes, collapse = ", "), "\\n")
} else {
  selected_celltypes <- LINEAGE_DEFINITIONS[[LINEAGE_TYPE]]
  cat("Selected lineage:", LINEAGE_TYPE, "\\n")
  cat("Cell types included:", paste(selected_celltypes, collapse = ", "), "\\n")
}

# Check which cell types are available
available_types <- unique(sobj@meta.data[[cell_type_col]])
missing_types <- setdiff(selected_celltypes, available_types)
present_types <- intersect(selected_celltypes, available_types)

if (length(missing_types) > 0) {
  cat("Warning: Cell types not found in dataset:", paste(missing_types, collapse = ", "), "\\n")
}

if (length(present_types) == 0) {
  stop("No selected cell types found in the dataset")
}

cat("Cell types found:", paste(present_types, collapse = ", "), "\\n")

# Subset to selected cell types
subset_condition <- paste0(cell_type_col, " %in% c('", paste(present_types, collapse = "', '"), "')")
sobj_subset <- subset(sobj, subset = eval(parse(text = subset_condition)))

cat("Subset contains", ncol(sobj_subset), "cells\\n")

# Set factor levels for consistent ordering
sobj_subset@meta.data[[cell_type_col]] <- factor(sobj_subset@meta.data[[cell_type_col]], 
                                                 levels = present_types)

# ============================================================================
# Quality Control and Sample Filtering
# ============================================================================

cat("Performing quality control...\\n")

# Check sample representation
sample_col <- "sample_id"
if (!"sample_id" %in% colnames(sobj_subset@meta.data)) {
  if ("id" %in% colnames(sobj_subset@meta.data)) {
    sample_col <- "id"
  } else if ("orig.ident" %in% colnames(sobj_subset@meta.data)) {
    sample_col <- "orig.ident"
  }
}

# Filter samples with sufficient cells
sample_counts <- table(sobj_subset@meta.data[[sample_col]])
valid_samples <- names(sample_counts)[sample_counts >= MIN_CELLS_PER_SAMPLE]

cat("Samples with >=", MIN_CELLS_PER_SAMPLE, "cells:", length(valid_samples), "/", length(sample_counts), "\\n")

if (length(valid_samples) == 0) {
  stop("No samples have sufficient cells for analysis")
}

# Filter to valid samples
if (length(valid_samples) < length(sample_counts)) {
  subset_condition <- paste0(sample_col, " %in% c('", paste(valid_samples, collapse = "', '"), "')")
  sobj_subset <- subset(sobj_subset, subset = eval(parse(text = subset_condition)))
  cat("After sample filtering:", ncol(sobj_subset), "cells\\n")
}

# Print cell type distribution
cat("\\nCell type distribution:\\n")
print(table(sobj_subset@meta.data[[cell_type_col]]))

cat("\\nSample distribution:\\n")
print(table(sobj_subset@meta.data[[sample_col]]))

# ============================================================================
# Data Preparation for Integration
# ============================================================================

cat("\\nPreparing data for integration...\\n")

# Extract count matrix
counts <- LayerData(sobj_subset, assay = "RNA", layer = "counts")

# Create new Seurat object
sobj_new <- CreateSeuratObject(counts = counts, meta.data = sobj_subset@meta.data)

# Split RNA assay by sample for integration
sobj_new[["RNA"]] <- split(sobj_new[["RNA"]], f = sobj_new@meta.data[[sample_col]])

cat("Split RNA assay by", length(unique(sobj_new@meta.data[[sample_col]])), "samples\\n")

# ============================================================================
# Standard Preprocessing
# ============================================================================

cat("Performing standard preprocessing...\\n")

# Normalization, variable features, scaling, and PCA
sobj_new <- NormalizeData(sobj_new, verbose = FALSE)
sobj_new <- FindVariableFeatures(sobj_new, verbose = FALSE)
sobj_new <- ScaleData(sobj_new, verbose = FALSE)
sobj_new <- RunPCA(sobj_new, verbose = FALSE)

cat("Found", length(VariableFeatures(sobj_new)), "variable features\\n")

# ============================================================================
# Integration Analysis
# ============================================================================

cat("Performing", INTEGRATION_METHOD, "integration...\\n")

# Construct reduction name
reduction_name <- paste('integrated', IDENTIFIER, sep = '.')

# Perform integration
sobj_new <- IntegrateLayers(
  object = sobj_new, 
  method = INTEGRATION_METHOD,
  new.reduction = reduction_name,
  verbose = FALSE
)

# Re-join layers after integration
sobj_new[["RNA"]] <- JoinLayers(sobj_new[["RNA"]])

cat("Integration completed with reduction:", reduction_name, "\\n")

# ============================================================================
# Clustering and Dimensionality Reduction
# ============================================================================

cat("Performing clustering and UMAP...\\n")

# Find neighbors and clusters
sobj_new <- FindNeighbors(sobj_new, reduction = reduction_name, dims = PCA_DIMS, verbose = FALSE)
sobj_new <- FindClusters(sobj_new, resolution = CLUSTERING_RESOLUTIONS, verbose = FALSE)

# UMAP embedding
sobj_new <- RunUMAP(sobj_new, dims = PCA_DIMS, reduction = reduction_name, verbose = FALSE)

# Report clustering results
n_clusters_default <- length(unique(sobj_new$seurat_clusters))
cat("Generated", n_clusters_default, "clusters at default resolution\\n")

# ============================================================================
# Visualization and Quality Assessment
# ============================================================================

cat("Creating visualizations...\\n")

# Original cell type UMAP
original_plot <- DimPlot(sobj_new, group.by = cell_type_col, label = TRUE, repel = TRUE, 
                        pt.size = 0.5) +
  ggtitle(paste("Original", LINEAGE_TYPE, "Cell Types")) +
  theme(legend.position = "bottom")

# New clusters UMAP
cluster_plot <- DimPlot(sobj_new, group.by = "seurat_clusters", label = TRUE, 
                       pt.size = 0.5) +
  ggtitle(paste("Subclustering Results (", IDENTIFIER, ")", sep = ""))

# Sample split plot for integration quality
sample_plot <- DimPlot(sobj_new, group.by = "seurat_clusters", split.by = sample_col, 
                      ncol = min(4, length(unique(sobj_new@meta.data[[sample_col]]))), 
                      pt.size = 0.3) +
  ggtitle(paste("Integration Quality Check -", IDENTIFIER))

# Multiple resolution comparison
if (length(CLUSTERING_RESOLUTIONS) >= 4) {
  # Select 4 representative resolutions
  res_indices <- round(seq(1, length(CLUSTERING_RESOLUTIONS), length.out = 4))
  selected_res <- CLUSTERING_RESOLUTIONS[res_indices]
  
  res_plots <- list()
  for (i in seq_along(selected_res)) {
    res_col <- paste0("RNA_snn_res.", selected_res[i])
    if (res_col %in% colnames(sobj_new@meta.data)) {
      res_plots[[i]] <- DimPlot(sobj_new, group.by = res_col, label = TRUE, pt.size = 0.3) +
        ggtitle(paste("Resolution", selected_res[i])) +
        theme(legend.position = "none")
    }
  }
  
  if (length(res_plots) > 0) {
    resolution_comparison <- wrap_plots(res_plots, ncol = 2)
  }
}

# Save main visualizations
main_plots <- original_plot / cluster_plot
ggsave(file.path(PLOTS_DIR, paste0("subclustering_", LINEAGE_TYPE, "_", IDENTIFIER, "_main.pdf")), 
       main_plots, width = 12, height = 12)

ggsave(file.path(PLOTS_DIR, paste0("subclustering_", LINEAGE_TYPE, "_", IDENTIFIER, "_quality.pdf")), 
       sample_plot, width = 16, height = 10)

if (exists("resolution_comparison")) {
  ggsave(file.path(PLOTS_DIR, paste0("subclustering_", LINEAGE_TYPE, "_", IDENTIFIER, "_resolutions.pdf")), 
         resolution_comparison, width = 12, height = 10)
}

# ============================================================================
# Cluster Characterization
# ============================================================================

cat("Characterizing subclusters...\\n")

# Find markers for new subclusters
Idents(sobj_new) <- "seurat_clusters"
subcluster_markers <- FindAllMarkers(sobj_new, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25, verbose = FALSE)

# Create cluster summary
cluster_summary <- sobj_new@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    n_samples = length(unique(get(sample_col))),
    dominant_celltype = names(sort(table(get(cell_type_col)), decreasing = TRUE))[1],
    celltype_purity = max(table(get(cell_type_col))) / n(),
    .groups = 'drop'
  )

# Add top markers to summary
if (nrow(subcluster_markers) > 0) {
  top_markers_per_cluster <- subcluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) %>%
    summarise(top_markers = paste(gene, collapse = ", "), .groups = 'drop')
  
  cluster_summary <- left_join(cluster_summary, top_markers_per_cluster, 
                              by = c("seurat_clusters" = "cluster"))
}

# ============================================================================
# Generate Output Files
# ============================================================================

cat("Generating output files...\\n")

# Define output file names
output_prefix <- paste0("subclustering_", LINEAGE_TYPE, "_", IDENTIFIER)
output_rds <- file.path(OUTPUT_DIR, paste0(output_prefix, ".rds"))
output_loupe <- file.path(OUTPUT_DIR, paste0(output_prefix, "_loupe"))

# Save Seurat object
saveRDS(sobj_new, output_rds)

# Save marker genes
if (nrow(subcluster_markers) > 0) {
  write.csv(subcluster_markers, file.path(OUTPUT_DIR, paste0(output_prefix, "_markers.csv")), 
            row.names = FALSE)
}

# Save cluster summary
write.csv(cluster_summary, file.path(OUTPUT_DIR, paste0(output_prefix, "_summary.csv")), 
          row.names = FALSE)

# Generate Loupe file
tryCatch({
  create_loupe_from_seurat(sobj_new, output_name = output_loupe, force = TRUE)
  cat("Loupe file created successfully\\n")
}, error = function(e) {
  cat("Warning: Could not create Loupe file:", e$message, "\\n")
})

# ============================================================================
# Generate Analysis Report
# ============================================================================

cat("Generating analysis report...\\n")

# Calculate summary statistics
n_original_types <- length(present_types)
n_subclusters <- length(unique(sobj_new$seurat_clusters))
n_total_cells <- ncol(sobj_new)
n_samples_used <- length(unique(sobj_new@meta.data[[sample_col]]))

# Generate comprehensive report
report_text <- paste0(
  "Lineage Subclustering Analysis Report\\n",
  "===================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Lineage: ", LINEAGE_TYPE, "\\n",
  "Integration Method: ", INTEGRATION_METHOD, "\\n\\n",
  "Data Summary:\\n",
  "  - Original cell types: ", n_original_types, "\\n",
  "  - Cell types included: ", paste(present_types, collapse = ", "), "\\n",
  "  - Total cells analyzed: ", n_total_cells, "\\n",
  "  - Samples included: ", n_samples_used, "\\n",
  "  - Subclusters identified: ", n_subclusters, "\\n\\n",
  "Analysis Parameters:\\n",
  "  - Integration method: ", INTEGRATION_METHOD, "\\n",
  "  - PCA dimensions: ", paste(range(PCA_DIMS), collapse = "-"), "\\n",
  "  - Clustering resolutions: ", paste(range(CLUSTERING_RESOLUTIONS), collapse = "-"), "\\n",
  "  - Min cells per sample: ", MIN_CELLS_PER_SAMPLE, "\\n\\n",
  "Results:\\n"
)

# Add cluster summary to report
for (i in 1:nrow(cluster_summary)) {
  cluster_info <- cluster_summary[i, ]
  report_text <- paste0(report_text,
    "  Subcluster ", cluster_info$seurat_clusters, ":\\n",
    "    - Cells: ", cluster_info$n_cells, "\\n",
    "    - Samples: ", cluster_info$n_samples, "\\n",
    "    - Dominant type: ", cluster_info$dominant_celltype, 
    " (", round(100 * cluster_info$celltype_purity, 1), "%)\\n"
  )
  
  if ("top_markers" %in% colnames(cluster_summary)) {
    report_text <- paste0(report_text,
      "    - Top markers: ", cluster_info$top_markers, "\\n"
    )
  }
  report_text <- paste0(report_text, "\\n")
}

report_text <- paste0(report_text,
  "Output Files:\\n",
  "  - Seurat object: ", basename(output_rds), "\\n",
  "  - Cluster summary: ", paste0(output_prefix, "_summary.csv"), "\\n"
)

if (nrow(subcluster_markers) > 0) {
  report_text <- paste0(report_text,
    "  - Marker genes: ", paste0(output_prefix, "_markers.csv"), "\\n"
  )
}

report_text <- paste0(report_text,
  "  - Main plots: ", paste0(output_prefix, "_main.pdf"), "\\n",
  "  - Quality plots: ", paste0(output_prefix, "_quality.pdf"), "\\n",
  "  - Loupe file: ", paste0(output_prefix, "_loupe"), "\\n\\n",
  "Recommendations:\\n",
  "1. Review subclusters for biological relevance\\n",
  "2. Validate marker genes with literature\\n",
  "3. Consider merging similar subclusters if appropriate\\n",
  "4. Use results for downstream trajectory or functional analysis\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, paste0(output_prefix, "_report.txt")))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, paste0(output_prefix, "_sessionInfo.txt")))

# ============================================================================
# Final Summary
# ============================================================================

cat("\\nSubclustering analysis completed successfully!\\n")
cat("Lineage analyzed:", LINEAGE_TYPE, "\\n")
cat("Integration method:", INTEGRATION_METHOD, "\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nSummary:\\n")
cat("- Original cell types:", n_original_types, "\\n")
cat("- Cells analyzed:", n_total_cells, "\\n")
cat("- Subclusters found:", n_subclusters, "\\n")
cat("- Samples included:", n_samples_used, "\\n")

if (nrow(subcluster_markers) > 0) {
  cat("- Marker genes identified:", nrow(subcluster_markers), "\\n")
}

cat("\\nUsage examples:\\n")
cat("Rscript lineage_subclustering.R 'Neural_lineage' 'HarmonyIntegration'\\n")
cat("Rscript lineage_subclustering.R 'Brain_and_spinal_cord,Neural_ectoderm' 'RPCAIntegration'\\n")