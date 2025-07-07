#!/usr/bin/env Rscript

# ============================================================================
# Marker Gene Analysis and Cell Type Annotation
# ============================================================================
# 
# This script identifies cluster-specific marker genes and assigns cell type
# annotations based on known marker gene expression patterns. It is designed
# for spatial transcriptomics data from primate embryos.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(writexl)
  library(ggplot2)
  library(patchwork)
})

# Parameters for marker gene analysis
MIN_PCT <- 0.25           # Minimum percentage of cells expressing the gene
LOGFC_THRESHOLD <- 0.25   # Minimum log fold change threshold
TEST_USE <- "wilcox"      # Statistical test for DE analysis
TOP_MARKERS <- 10         # Number of top markers to display per cluster

# Input/Output paths
INPUT_SEURAT <- "../output/sobj_initial_clustering.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Known Marker Genes for Primate Embryo Cell Types
# ============================================================================

# Define known marker genes for major lineages based on literature
MARKER_GENES <- list(
  # Neural/Ectodermal lineages
  "Brain_and_spinal_cord" = c("SOX2", "SFRP2", "PAX6", "OLIG2", "FOXG1"),
  "Floor_plate" = c("FOXB1", "NKX6-1", "NKX6-2", "SLIT1", "SLIT2"),
  
  # Mesodermal lineages  
  "Paraxial_mesoderm" = c("MEOX1", "TCF15", "EGR2", "CRABP1", "ALDH1A2"),
  "Caudal_mesoderm" = c("HES7", "CDX4", "TBX6", "EVX1", "TBXT"),
  "Lateral_plate_mesoderm" = c("LIX1", "MEIS2", "HAND1", "HAND2"),
  "Heart" = c("MYH6", "TTN", "NKX2-5", "GATA4", "ACTC1", "NPPA"),
  "Cardiac_progenitors" = c("ISL1", "IRX3", "SALL1", "HES1", "MYCN"),
  "Endothelial_cells" = c("SOX7", "PECAM1", "CDH5", "NFATC1", "KLF2"),
  
  # Endodermal lineages
  "Gut_tube" = c("FOXA1", "FOXA2", "EPHB1", "SOX17"),
  "Foregut" = c("CALB1", "HHEX", "TTR", "HNF4A", "ONECUT1"),
  "Midgut" = c("PRSS2", "FOXJ1", "CDH6", "PBX1"),
  "Hindgut" = c("CDX4", "CDX2", "CDX1", "FOXA1", "HNF1B"),
  
  # Axial structures
  "Notochord" = c("NOTO", "FOXA2", "CHRD", "WNT5B", "NKX1-2", "FZD10"),
  "Prechordal_plate" = c("OTX2", "GSC", "TBXT", "IRX2", "DKK1", "NOG", "FST"),
  
  # Extraembryonic tissues
  "Allantois" = c("PENK", "PITX1", "PITX2", "HOXA10", "HOXA11", "HOXA13"),
  "Amnion" = c("GABRP", "POSTN"),
  "Yolk_sac_endoderm" = c("APOM", "SERPINA1", "AGT", "GATA1"),
  "Yolk_sac_mesoderm" = c("CD44", "AQP1"),
  "Blood_progenitor" = c("GP9", "THBS1", "PLEK", "KLF1"),
  
  # Signaling and development
  "PGC_markers" = c("NANOG", "NANOS3", "TCL1B", "TFAP2C"),
  "Chemokine_signaling" = c("CXCR4", "CXCR7", "CXCL12")
)

# ============================================================================
# Load Data and Prepare for Analysis
# ============================================================================

cat("Loading Seurat object...\n")

# Load the clustered Seurat object
if (!file.exists(INPUT_SEURAT)) {
  stop("Clustered Seurat object not found. Please run initial_clustering.R first.")
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded object with", ncol(sobj), "cells and", nrow(sobj), "features\n")
cat("Number of clusters:", length(unique(Idents(sobj))), "\n")

# ============================================================================
# Find Cluster Markers
# ============================================================================

cat("Finding cluster-specific marker genes...\n")

# Find markers for all clusters
all_markers <- FindAllMarkers(
  sobj,
  only.pos = TRUE,
  min.pct = MIN_PCT,
  logfc.threshold = LOGFC_THRESHOLD,
  test.use = TEST_USE,
  verbose = FALSE
)

# Filter and rank markers
all_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

cat("Found", nrow(all_markers), "significant marker genes\n")

# Get top markers per cluster
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = TOP_MARKERS, wt = avg_log2FC) %>%
  ungroup()

# Save all markers
write.csv(all_markers, file.path(OUTPUT_DIR, "all_cluster_markers.csv"), 
          row.names = FALSE)
write.csv(top_markers, file.path(OUTPUT_DIR, "top_cluster_markers.csv"), 
          row.names = FALSE)

# ============================================================================
# Visualize Top Markers
# ============================================================================

cat("Creating marker gene visualizations...\n")

# Heatmap of top markers
if (nrow(top_markers) > 0) {
  top_marker_genes <- unique(top_markers$gene)
  
  # Limit to reasonable number for visualization
  if (length(top_marker_genes) > 100) {
    top_marker_genes <- head(top_marker_genes, 100)
  }
  
  # Create heatmap
  marker_heatmap <- DoHeatmap(sobj, features = top_marker_genes, size = 3) +
    theme(axis.text.y = element_text(size = 6))
  
  ggsave(file.path(PLOTS_DIR, "marker_genes_heatmap.pdf"), 
         marker_heatmap, width = 12, height = 20, limitsize = FALSE)
  
  # Dot plot of selected markers
  if (length(top_marker_genes) <= 50) {
    marker_dotplot <- DotPlot(sobj, features = top_marker_genes) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(PLOTS_DIR, "marker_genes_dotplot.pdf"), 
           marker_dotplot, width = 15, height = 8)
  }
}

# ============================================================================
# Known Marker Gene Analysis
# ============================================================================

cat("Analyzing expression of known marker genes...\n")

# Flatten the marker gene list
all_known_markers <- unique(unlist(MARKER_GENES))

# Check which markers are present in the dataset
present_markers <- all_known_markers[all_known_markers %in% rownames(sobj)]
missing_markers <- all_known_markers[!all_known_markers %in% rownames(sobj)]

cat("Known markers present in dataset:", length(present_markers), "/", 
    length(all_known_markers), "\n")

if (length(missing_markers) > 0) {
  cat("Missing markers:", paste(head(missing_markers, 10), collapse = ", "), 
      ifelse(length(missing_markers) > 10, "...", ""), "\n")
}

# Create feature plots for key markers
if (length(present_markers) > 0) {
  # Select representative markers for visualization
  key_markers <- c("SOX2", "MEOX1", "NOTO", "MYH6", "FOXA1", "PECAM1", 
                   "CDX4", "APOM", "PENK")
  key_markers <- key_markers[key_markers %in% present_markers]
  
  if (length(key_markers) > 0) {
    # UMAP feature plots
    feature_plots <- list()
    for (i in seq_along(key_markers)) {
      feature_plots[[i]] <- FeaturePlot(sobj, features = key_markers[i], 
                                        pt.size = 0.3) +
        theme(axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank())
    }
    
    # Combine plots
    combined_features <- wrap_plots(feature_plots, ncol = 3)
    ggsave(file.path(PLOTS_DIR, "key_markers_featureplot.pdf"), 
           combined_features, width = 15, height = ceiling(length(key_markers)/3) * 4)
    
    # Violin plots
    violin_plot <- VlnPlot(sobj, features = key_markers, pt.size = 0.1, ncol = 3)
    ggsave(file.path(PLOTS_DIR, "key_markers_violin.pdf"), 
           violin_plot, width = 15, height = ceiling(length(key_markers)/3) * 4)
  }
}

# ============================================================================
# Automated Cell Type Annotation
# ============================================================================

cat("Performing automated cell type annotation...\n")

# Calculate average expression of marker gene sets for each cluster
cluster_scores <- data.frame(cluster = levels(Idents(sobj)))

for (cell_type in names(MARKER_GENES)) {
  markers <- MARKER_GENES[[cell_type]]
  markers_present <- markers[markers %in% rownames(sobj)]
  
  if (length(markers_present) > 0) {
    # Calculate average expression per cluster
    avg_exp <- AverageExpression(sobj, features = markers_present, 
                                group.by = "seurat_clusters")$RNA
    
    if (length(markers_present) == 1) {
      cluster_scores[[cell_type]] <- as.numeric(avg_exp)
    } else {
      cluster_scores[[cell_type]] <- colMeans(avg_exp, na.rm = TRUE)
    }
  } else {
    cluster_scores[[cell_type]] <- 0
  }
}

# Find best matching cell type for each cluster
cluster_annotations <- cluster_scores %>%
  column_to_rownames("cluster") %>%
  apply(1, function(x) {
    if (all(x == 0)) return("Unknown")
    names(x)[which.max(x)]
  })

# Create annotation data frame
annotation_df <- data.frame(
  cluster = names(cluster_annotations),
  predicted_celltype = cluster_annotations,
  stringsAsFactors = FALSE
)

# Add confidence scores (max score / second max score ratio)
confidence_scores <- cluster_scores %>%
  column_to_rownames("cluster") %>%
  apply(1, function(x) {
    if (all(x == 0)) return(0)
    sorted_x <- sort(x, decreasing = TRUE)
    if (length(sorted_x) == 1 || sorted_x[2] == 0) return(Inf)
    sorted_x[1] / sorted_x[2]
  })

annotation_df$confidence_score <- confidence_scores[annotation_df$cluster]

# Save annotation results
write.csv(annotation_df, file.path(OUTPUT_DIR, "automated_annotations.csv"), 
          row.names = FALSE)
write.csv(cluster_scores, file.path(OUTPUT_DIR, "cluster_marker_scores.csv"), 
          row.names = TRUE)

# ============================================================================
# Manual Annotation Guidelines
# ============================================================================

cat("Generating manual annotation guidelines...\n")

# Create a summary table for manual review
manual_review <- all_markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    top_3_markers = paste(head(gene, 3), collapse = ", "),
    avg_logFC_range = paste(round(range(avg_log2FC), 2), collapse = " - "),
    .groups = 'drop'
  ) %>%
  left_join(annotation_df, by = "cluster")

write.csv(manual_review, file.path(OUTPUT_DIR, "manual_annotation_guide.csv"), 
          row.names = FALSE)

# ============================================================================
# Apply Preliminary Annotations
# ============================================================================

cat("Applying preliminary annotations to Seurat object...\n")

# Add annotations to Seurat object
sobj$preliminary_celltype <- annotation_df$predicted_celltype[
  match(as.character(Idents(sobj)), annotation_df$cluster)
]

# Create annotated plots
annotated_umap <- DimPlot(sobj, group.by = "preliminary_celltype", 
                         label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("Preliminary Cell Type Annotations") +
  theme(legend.position = "bottom")

ggsave(file.path(PLOTS_DIR, "preliminary_annotations_umap.pdf"), 
       annotated_umap, width = 12, height = 10)

# Compare with cluster plot
cluster_umap <- DimPlot(sobj, label = TRUE, pt.size = 0.3) +
  ggtitle("Original Clusters")

comparison_plot <- cluster_umap | annotated_umap
ggsave(file.path(PLOTS_DIR, "cluster_vs_annotation_comparison.pdf"), 
       comparison_plot, width = 20, height = 8)

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating comprehensive annotation report...\n")

# Summary statistics
n_clusters <- length(unique(Idents(sobj)))
n_annotated <- sum(annotation_df$predicted_celltype != "Unknown")
n_high_confidence <- sum(annotation_df$confidence_score > 2, na.rm = TRUE)

# Create detailed report
report_text <- paste0(
  "Marker Gene Analysis and Cell Type Annotation Report\n",
  "===================================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Dataset: ", ncol(sobj), " cells, ", nrow(sobj), " features\n\n",
  "Cluster Analysis:\n",
  "  - Total clusters: ", n_clusters, "\n",
  "  - Clusters with annotations: ", n_annotated, "\n",
  "  - High confidence annotations: ", n_high_confidence, "\n\n",
  "Marker Gene Analysis:\n",
  "  - Total significant markers found: ", nrow(all_markers), "\n",
  "  - Known markers in dataset: ", length(present_markers), "/", 
      length(all_known_markers), "\n",
  "  - Statistical test used: ", TEST_USE, "\n",
  "  - Min pct threshold: ", MIN_PCT, "\n",
  "  - LogFC threshold: ", LOGFC_THRESHOLD, "\n\n",
  "Output Files:\n",
  "  - all_cluster_markers.csv: Complete marker gene list\n",
  "  - top_cluster_markers.csv: Top markers per cluster\n",
  "  - automated_annotations.csv: Predicted cell types\n",
  "  - cluster_marker_scores.csv: Expression scores for known markers\n",
  "  - manual_annotation_guide.csv: Guide for manual review\n\n",
  "Visualization Files:\n",
  "  - marker_genes_heatmap.pdf: Heatmap of top markers\n",
  "  - key_markers_featureplot.pdf: Expression patterns of key markers\n",
  "  - preliminary_annotations_umap.pdf: Annotated UMAP\n",
  "  - cluster_vs_annotation_comparison.pdf: Cluster vs annotation comparison\n\n",
  "Next Steps:\n",
  "1. Review automated annotations using manual_annotation_guide.csv\n",
  "2. Validate annotations using spatial information and literature\n",
  "3. Refine annotations based on biological knowledge\n",
  "4. Consider subclustering for complex lineages\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "annotation_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving annotated Seurat object...\n")

# Save the object with preliminary annotations
saveRDS(sobj, file.path(OUTPUT_DIR, "sobj_with_preliminary_annotations.rds"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_marker_annotation.txt"))

cat("\nMarker gene analysis and annotation completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")
cat("\nPlease review the automated annotations and refine as needed based on:\n")
cat("1. Spatial information\n")
cat("2. Literature knowledge\n") 
cat("3. Manual inspection of marker gene expression\n")

# Print summary table
cat("\nPreliminary Annotation Summary:\n")
print(table(sobj$preliminary_celltype))