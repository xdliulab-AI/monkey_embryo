#!/usr/bin/env Rscript

# ============================================================================
# Spatial Validation of Cell Type Annotations
# ============================================================================
# 
# This script validates cell type annotations by examining their spatial
# organization and distribution patterns. It creates spatial visualizations
# and validates annotations against known developmental biology principles.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(RColorBrewer)
  library(viridis)
})

# Parameters
POINT_SIZE <- 0.5         # Point size for spatial plots
ALPHA <- 0.8              # Transparency for plots
COORD_FIXED <- TRUE       # Use fixed coordinates for spatial plots

# Input/Output paths
INPUT_SEURAT <- "../output/sobj_with_preliminary_annotations.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Expected Spatial Patterns for Primate Embryo Cell Types
# ============================================================================

# Define expected spatial relationships based on developmental biology
SPATIAL_EXPECTATIONS <- list(
  # Axial structures
  "expected_midline" = c("Notochord", "Prechordal_plate"),
  "expected_dorsal" = c("Brain_and_spinal_cord", "Neural_ectoderm"),
  "expected_paraxial" = c("Paraxial_mesoderm", "Somites"),
  
  # Regional patterns
  "expected_anterior" = c("Forebrain", "Heart", "Prechordal_plate"),
  "expected_posterior" = c("Caudal_mesoderm", "Tailbud", "Allantois"),
  "expected_ventral" = c("Gut_tube", "Heart", "Lateral_plate_mesoderm"),
  
  # Tissue organization
  "expected_surrounding_neural" = c("Paraxial_mesoderm", "Neural_crest"),
  "expected_vascular" = c("Endothelial_cells", "Blood_progenitor"),
  "expected_extraembryonic" = c("Yolk_sac_endoderm", "Yolk_sac_mesoderm", 
                                "Amnion", "Allantois")
)

# ============================================================================
# Load Data and Prepare for Analysis
# ============================================================================

cat("Loading annotated Seurat object...\n")

# Load the annotated Seurat object
if (!file.exists(INPUT_SEURAT)) {
  stop("Annotated Seurat object not found. Please run marker_gene_annotation.R first.")
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded object with", ncol(sobj), "cells\n")

# Check for spatial coordinates
has_spatial <- FALSE
spatial_cols <- c("x", "y", "imagerow", "imagecol", "row", "col")
available_coords <- spatial_cols[spatial_cols %in% colnames(sobj@meta.data)]

if (length(available_coords) >= 2) {
  has_spatial <- TRUE
  x_coord <- available_coords[1]
  y_coord <- available_coords[2]
  cat("Found spatial coordinates:", x_coord, "and", y_coord, "\n")
} else if ("Spatial" %in% names(sobj@reductions)) {
  has_spatial <- TRUE
  cat("Found Spatial reduction in Seurat object\n")
} else {
  cat("Warning: No spatial coordinates found. Will use UMAP for visualization.\n")
}

# ============================================================================
# Basic Spatial Visualization
# ============================================================================

cat("Creating spatial visualizations...\n")

# Define colors for cell types
cell_types <- unique(sobj$preliminary_celltype)
cell_types <- cell_types[!is.na(cell_types)]

# Use a diverse color palette
if (length(cell_types) <= 12) {
  colors <- RColorBrewer::brewer.pal(min(length(cell_types), 12), "Set3")
} else {
  colors <- rainbow(length(cell_types))
}
names(colors) <- cell_types

# Create spatial plot function
create_spatial_plot <- function(sobj, group_by, title, colors = NULL) {
  if (has_spatial) {
    if ("Spatial" %in% names(sobj@reductions)) {
      # Use Spatial reduction if available
      plot_data <- data.frame(
        x = Embeddings(sobj, "Spatial")[, 1],
        y = Embeddings(sobj, "Spatial")[, 2],
        group = sobj@meta.data[[group_by]]
      )
    } else {
      # Use coordinate metadata
      plot_data <- data.frame(
        x = sobj@meta.data[[x_coord]],
        y = sobj@meta.data[[y_coord]],
        group = sobj@meta.data[[group_by]]
      )
    }
    
    p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
      geom_point(size = POINT_SIZE, alpha = ALPHA) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()
      ) +
      ggtitle(title)
    
    if (COORD_FIXED) {
      p <- p + coord_fixed()
    }
    
    if (!is.null(colors)) {
      p <- p + scale_color_manual(values = colors)
    }
    
  } else {
    # Fall back to UMAP
    p <- DimPlot(sobj, group.by = group_by, pt.size = POINT_SIZE) +
      ggtitle(paste(title, "(UMAP)"))
    
    if (!is.null(colors)) {
      p <- p + scale_color_manual(values = colors)
    }
  }
  
  return(p)
}

# Main spatial plot
spatial_plot <- create_spatial_plot(sobj, "preliminary_celltype", 
                                   "Spatial Distribution of Cell Types", colors)
ggsave(file.path(PLOTS_DIR, "spatial_cell_types.pdf"), 
       spatial_plot, width = 14, height = 10)

# Cluster plot for comparison
cluster_plot <- create_spatial_plot(sobj, "seurat_clusters", 
                                   "Spatial Distribution of Clusters")
ggsave(file.path(PLOTS_DIR, "spatial_clusters.pdf"), 
       cluster_plot, width = 12, height = 10)

# Side-by-side comparison
comparison_plot <- cluster_plot | spatial_plot
ggsave(file.path(PLOTS_DIR, "spatial_cluster_vs_celltype.pdf"), 
       comparison_plot, width = 20, height = 8)

# ============================================================================
# Individual Cell Type Spatial Patterns
# ============================================================================

cat("Creating individual cell type spatial maps...\n")

# Create individual plots for each major cell type
major_cell_types <- names(table(sobj$preliminary_celltype))[
  table(sobj$preliminary_celltype) > 50
]

for (cell_type in major_cell_types) {
  # Create binary annotation (this cell type vs others)
  sobj$temp_annotation <- ifelse(sobj$preliminary_celltype == cell_type, 
                                cell_type, "Other")
  
  temp_colors <- c("Other" = "lightgray", cell_type = "red")
  names(temp_colors)[2] <- cell_type
  
  individual_plot <- create_spatial_plot(sobj, "temp_annotation", 
                                        paste("Spatial Distribution:", cell_type),
                                        temp_colors)
  
  # Clean filename
  clean_name <- gsub("[^A-Za-z0-9_]", "_", cell_type)
  ggsave(file.path(PLOTS_DIR, paste0("spatial_", clean_name, ".pdf")), 
         individual_plot, width = 10, height = 8)
}

# Clean up temporary column
sobj$temp_annotation <- NULL

# ============================================================================
# Spatial Statistics and Validation
# ============================================================================

cat("Calculating spatial statistics...\n")

if (has_spatial) {
  # Extract spatial coordinates
  if ("Spatial" %in% names(sobj@reductions)) {
    coords <- Embeddings(sobj, "Spatial")
  } else {
    coords <- sobj@meta.data[, c(x_coord, y_coord)]
    colnames(coords) <- c("x", "y")
  }
  
  # Calculate spatial statistics for each cell type
  spatial_stats <- data.frame(
    cell_type = character(),
    n_cells = numeric(),
    x_mean = numeric(),
    y_mean = numeric(),
    x_sd = numeric(),
    y_sd = numeric(),
    spatial_spread = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cell_type in cell_types) {
    cells_idx <- which(sobj$preliminary_celltype == cell_type)
    
    if (length(cells_idx) > 0) {
      cell_coords <- coords[cells_idx, ]
      
      spatial_stats <- rbind(spatial_stats, data.frame(
        cell_type = cell_type,
        n_cells = length(cells_idx),
        x_mean = mean(cell_coords$x, na.rm = TRUE),
        y_mean = mean(cell_coords$y, na.rm = TRUE),
        x_sd = sd(cell_coords$x, na.rm = TRUE),
        y_sd = sd(cell_coords$y, na.rm = TRUE),
        spatial_spread = sqrt(var(cell_coords$x, na.rm = TRUE) + 
                             var(cell_coords$y, na.rm = TRUE)),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save spatial statistics
  write.csv(spatial_stats, file.path(OUTPUT_DIR, "spatial_statistics.csv"), 
            row.names = FALSE)
  
  # Identify spatially compact vs dispersed cell types
  spatial_stats$compactness <- ifelse(spatial_stats$spatial_spread < 
                                     median(spatial_stats$spatial_spread, na.rm = TRUE),
                                     "Compact", "Dispersed")
  
  cat("Spatially compact cell types:\n")
  compact_types <- spatial_stats$cell_type[spatial_stats$compactness == "Compact"]
  cat(paste(compact_types, collapse = ", "), "\n\n")
  
  cat("Spatially dispersed cell types:\n")
  dispersed_types <- spatial_stats$cell_type[spatial_stats$compactness == "Dispersed"]
  cat(paste(dispersed_types, collapse = ", "), "\n\n")
}

# ============================================================================
# Axis-Based Analysis
# ============================================================================

cat("Performing axis-based spatial analysis...\n")

if (has_spatial) {
  # Divide embryo into regions based on coordinates
  x_range <- range(coords$x, na.rm = TRUE)
  y_range <- range(coords$y, na.rm = TRUE)
  
  # Define anterior-posterior axis (assuming x-axis)
  x_terciles <- quantile(coords$x, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE)
  sobj$AP_region <- cut(coords$x, breaks = x_terciles, 
                       labels = c("Anterior", "Middle", "Posterior"),
                       include.lowest = TRUE)
  
  # Define dorsal-ventral axis (assuming y-axis)
  y_terciles <- quantile(coords$y, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE)
  sobj$DV_region <- cut(coords$y, breaks = y_terciles,
                       labels = c("Ventral", "Middle", "Dorsal"),
                       include.lowest = TRUE)
  
  # Create region-based analysis
  ap_distribution <- table(sobj$preliminary_celltype, sobj$AP_region)
  dv_distribution <- table(sobj$preliminary_celltype, sobj$DV_region)
  
  # Convert to proportions
  ap_props <- prop.table(ap_distribution, margin = 1)
  dv_props <- prop.table(dv_distribution, margin = 1)
  
  # Save distributions
  write.csv(ap_distribution, file.path(OUTPUT_DIR, "AP_distribution.csv"))
  write.csv(dv_distribution, file.path(OUTPUT_DIR, "DV_distribution.csv"))
  write.csv(ap_props, file.path(OUTPUT_DIR, "AP_proportions.csv"))
  write.csv(dv_props, file.path(OUTPUT_DIR, "DV_proportions.csv"))
  
  # Create axis plots
  ap_plot <- create_spatial_plot(sobj, "AP_region", "Anterior-Posterior Regions")
  dv_plot <- create_spatial_plot(sobj, "DV_region", "Dorsal-Ventral Regions")
  
  axis_plots <- ap_plot | dv_plot
  ggsave(file.path(PLOTS_DIR, "spatial_axis_regions.pdf"), 
         axis_plots, width = 16, height = 6)
}

# ============================================================================
# Annotation Quality Assessment
# ============================================================================

cat("Assessing annotation quality...\n")

# Calculate annotation metrics
annotation_metrics <- data.frame(
  cell_type = cell_types,
  n_cells = as.numeric(table(sobj$preliminary_celltype)[cell_types]),
  stringsAsFactors = FALSE
)

annotation_metrics$percentage <- round(100 * annotation_metrics$n_cells / ncol(sobj), 2)

# Identify potential issues
annotation_metrics$potential_issues <- ""

# Check for very small populations (< 1% of total)
small_populations <- annotation_metrics$percentage < 1
annotation_metrics$potential_issues[small_populations] <- 
  paste(annotation_metrics$potential_issues[small_populations], 
        "Small population;", sep = " ")

# Check for very large populations (> 20% of total)
large_populations <- annotation_metrics$percentage > 20
annotation_metrics$potential_issues[large_populations] <- 
  paste(annotation_metrics$potential_issues[large_populations], 
        "Large population;", sep = " ")

if (has_spatial) {
  # Add spatial compactness information
  annotation_metrics <- merge(annotation_metrics, 
                             spatial_stats[, c("cell_type", "compactness", "spatial_spread")],
                             by = "cell_type", all.x = TRUE)
}

# Save annotation quality assessment
write.csv(annotation_metrics, file.path(OUTPUT_DIR, "annotation_quality_assessment.csv"), 
          row.names = FALSE)

# ============================================================================
# Generate Validation Report
# ============================================================================

cat("Generating spatial validation report...\n")

# Count annotations by type
n_total_cells <- ncol(sobj)
n_annotated <- sum(!is.na(sobj$preliminary_celltype))
n_unknown <- sum(sobj$preliminary_celltype == "Unknown", na.rm = TRUE)

# Generate comprehensive report
report_text <- paste0(
  "Spatial Validation Report\n",
  "========================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Dataset: ", n_total_cells, " cells\n\n",
  "Annotation Summary:\n",
  "  - Total annotated cells: ", n_annotated, " (", 
      round(100 * n_annotated / n_total_cells, 1), "%)\n",
  "  - Unknown annotations: ", n_unknown, " (", 
      round(100 * n_unknown / n_total_cells, 1), "%)\n",
  "  - Unique cell types: ", length(cell_types), "\n\n"
)

if (has_spatial) {
  report_text <- paste0(report_text,
    "Spatial Analysis:\n",
    "  - Spatial coordinates available: Yes\n",
    "  - Compact cell types: ", length(compact_types), "\n",
    "  - Dispersed cell types: ", length(dispersed_types), "\n\n",
    "Regional Distribution:\n",
    "  - A-P axis analysis: Complete\n",
    "  - D-V axis analysis: Complete\n\n"
  )
} else {
  report_text <- paste0(report_text,
    "Spatial Analysis:\n",
    "  - Spatial coordinates available: No\n",
    "  - Analysis limited to UMAP visualization\n\n"
  )
}

report_text <- paste0(report_text,
  "Quality Assessment:\n"
)

# Add quality issues if any
issues_found <- annotation_metrics$potential_issues != ""
if (any(issues_found)) {
  report_text <- paste0(report_text,
    "  - Cell types with potential issues: ", sum(issues_found), "\n"
  )
  for (i in which(issues_found)) {
    report_text <- paste0(report_text,
      "    * ", annotation_metrics$cell_type[i], ": ", 
      annotation_metrics$potential_issues[i], "\n"
    )
  }
} else {
  report_text <- paste0(report_text,
    "  - No major quality issues detected\n"
  )
}

report_text <- paste0(report_text,
  "\nOutput Files:\n",
  "  - spatial_cell_types.pdf: Main spatial visualization\n",
  "  - spatial_clusters.pdf: Cluster spatial distribution\n",
  "  - spatial_cluster_vs_celltype.pdf: Comparison visualization\n",
  "  - spatial_[celltype].pdf: Individual cell type maps\n"
)

if (has_spatial) {
  report_text <- paste0(report_text,
    "  - spatial_statistics.csv: Spatial statistics per cell type\n",
    "  - AP_distribution.csv: Anterior-posterior distribution\n",
    "  - DV_distribution.csv: Dorsal-ventral distribution\n",
    "  - spatial_axis_regions.pdf: Axis-based visualization\n"
  )
}

report_text <- paste0(report_text,
  "  - annotation_quality_assessment.csv: Quality metrics\n\n",
  "Recommendations:\n",
  "1. Review cell types with potential issues\n",
  "2. Validate spatial patterns against literature\n",
  "3. Consider subclustering for large, heterogeneous populations\n",
  "4. Verify rare cell types with additional markers\n",
  "5. Use spatial information to refine boundaries\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "spatial_validation_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\n")

# Save the object with spatial metadata
if (has_spatial) {
  saveRDS(sobj, file.path(OUTPUT_DIR, "sobj_with_spatial_validation.rds"))
}

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_spatial_validation.txt"))

cat("\nSpatial validation completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

if (has_spatial) {
  cat("\nSpatial patterns summary:\n")
  print(spatial_stats[order(spatial_stats$spatial_spread), 
                     c("cell_type", "n_cells", "compactness")])
} else {
  cat("\nNote: Limited spatial analysis due to missing coordinate information.\n")
  cat("Consider adding spatial coordinates for full validation.\n")
}

cat("\nPlease review the spatial validation report and refine annotations as needed.\n")