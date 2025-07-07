#!/usr/bin/env Rscript

# ============================================================================
# Cross-Dataset Correlation and Integration Analysis
# ============================================================================
# 
# This script performs comprehensive cross-dataset correlation analysis to
# compare gene expression patterns between different datasets, developmental
# stages, or experimental conditions. It includes correlation heatmaps,
# integration quality assessment, and comparative analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(corrplot)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
  library(writexl)
})

# Parameters
JOB_ID <- "correlation_analysis"
DATASET_1_NAME <- "CS9"           # Name for first dataset
DATASET_2_NAME <- "CS10"          # Name for second dataset
CORRELATION_METHOD <- "spearman"   # "spearman" or "pearson"
TOP_GENES_OPTIONS <- c(25, 50, 100, 200, 500)  # Different numbers of top genes to test
CELL_TYPE_COL <- "celltype"       # Column name for cell types
CONDITION_COL <- "CS"             # Column name for conditions/datasets
MIN_CELLS_PER_TYPE <- 10          # Minimum cells per cell type

# Input/Output paths
INPUT_SEURAT <- "../data/annotated_seurat_object.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"
DATA_DIR <- "../data"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  INPUT_SEURAT <- args[1]
}
if (length(args) > 1) {
  DATASET_1_NAME <- args[2]
}
if (length(args) > 2) {
  DATASET_2_NAME <- args[3]
}
if (length(args) > 3) {
  JOB_ID <- args[4]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to prepare expression matrix for correlation analysis
prepare_expression_matrix <- function(sobj, condition_name, cell_type_col = CELL_TYPE_COL) {
  
  # Create combined identifier
  sobj$condition_celltype <- paste0(condition_name, "_", sobj@meta.data[[cell_type_col]])
  
  # Filter out cell types with too few cells
  cell_type_counts <- table(sobj@meta.data[[cell_type_col]])
  valid_cell_types <- names(cell_type_counts)[cell_type_counts >= MIN_CELLS_PER_TYPE]
  
  if (length(valid_cell_types) == 0) {
    stop("No cell types with sufficient cells found")
  }
  
  sobj_filtered <- subset(sobj, cells = which(sobj@meta.data[[cell_type_col]] %in% valid_cell_types))
  sobj_filtered$condition_celltype <- paste0(condition_name, "_", sobj_filtered@meta.data[[cell_type_col]])
  
  # Aggregate expression by condition_celltype
  avgexp <- AggregateExpression(
    sobj_filtered,
    return.seurat = TRUE,
    assays = "RNA",
    group.by = "condition_celltype"
  )
  
  # Extract expression matrix
  exp_matrix <- as.matrix(LayerData(avgexp, assay = "RNA", layer = "data"))
  
  return(exp_matrix)
}

# Function to load and process DEG data
load_deg_data <- function(deg_file) {
  if (file.exists(deg_file)) {
    deg_data <- readRDS(deg_file)
    cat("Loaded DEG data with", nrow(deg_data), "genes\n")
    return(deg_data)
  } else {
    cat("DEG file not found:", deg_file, "\n")
    return(NULL)
  }
}

# Function to select top genes by different criteria
select_top_genes <- function(deg_data, top_n, method = "avg_log2FC") {
  if (is.null(deg_data) || nrow(deg_data) == 0) {
    return(character(0))
  }
  
  if (method == "avg_log2FC") {
    top_genes <- deg_data %>%
      group_by(cluster) %>%
      arrange(desc(abs(avg_log2FC))) %>%
      slice_head(n = top_n) %>%
      pull(gene)
  } else if (method == "p_val_adj") {
    top_genes <- deg_data %>%
      group_by(cluster) %>%
      arrange(p_val_adj) %>%
      slice_head(n = top_n) %>%
      pull(gene)
  }
  
  return(unique(top_genes))
}

# Function to calculate correlation matrix
calculate_correlation_matrix <- function(mat1, mat2, common_genes, method = CORRELATION_METHOD) {
  
  if (length(common_genes) == 0) {
    warning("No common genes found for correlation analysis")
    return(NULL)
  }
  
  # Filter matrices by common genes
  mat1_filtered <- mat1[common_genes, , drop = FALSE]
  mat2_filtered <- mat2[common_genes, , drop = FALSE]
  
  # Combine matrices
  combined_matrix <- cbind(mat1_filtered, mat2_filtered)
  
  # Calculate correlation
  correlation_matrix <- cor(combined_matrix, method = method, use = "complete.obs")
  
  return(correlation_matrix)
}

# Function to create correlation heatmap
create_correlation_heatmap <- function(correlation_matrix, title, output_file, 
                                     dataset1_name = DATASET_1_NAME, 
                                     dataset2_name = DATASET_2_NAME) {
  
  if (is.null(correlation_matrix)) {
    cat("No correlation matrix to plot\n")
    return(NULL)
  }
  
  # Define color palette
  color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  breaks <- seq(-1, 1, length.out = 101)
  
  # Create annotation for datasets
  col_annotation <- data.frame(
    Dataset = ifelse(grepl(dataset1_name, colnames(correlation_matrix)), 
                    dataset1_name, dataset2_name)
  )
  rownames(col_annotation) <- colnames(correlation_matrix)
  
  # Define annotation colors
  ann_colors <- list(
    Dataset = c(setNames("red", dataset1_name), setNames("blue", dataset2_name))
  )
  
  # Create heatmap
  pdf(output_file, width = 12, height = 12)
  
  pheatmap(
    correlation_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = color_palette,
    breaks = breaks,
    main = title,
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 6,
    fontsize_row = 8,
    fontsize_col = 8,
    legend_breaks = seq(-1, 1, by = 0.5),
    legend_labels = seq(-1, 1, by = 0.5),
    annotation_col = col_annotation,
    annotation_colors = ann_colors,
    angle_col = 45
  )
  
  dev.off()
  
  cat("Correlation heatmap saved:", output_file, "\n")
  return(correlation_matrix)
}

# Function to calculate correlation statistics
calculate_correlation_stats <- function(correlation_matrix, dataset1_name, dataset2_name) {
  
  if (is.null(correlation_matrix)) {
    return(NULL)
  }
  
  # Identify cross-dataset correlations
  dataset1_cols <- grep(dataset1_name, colnames(correlation_matrix))
  dataset2_cols <- grep(dataset2_name, colnames(correlation_matrix))
  
  if (length(dataset1_cols) == 0 || length(dataset2_cols) == 0) {
    warning("Could not identify dataset columns for correlation statistics")
    return(NULL)
  }
  
  # Extract cross-dataset correlation submatrix
  cross_corr <- correlation_matrix[dataset1_cols, dataset2_cols, drop = FALSE]
  
  # Calculate statistics
  stats <- list(
    mean_correlation = mean(cross_corr, na.rm = TRUE),
    median_correlation = median(cross_corr, na.rm = TRUE),
    min_correlation = min(cross_corr, na.rm = TRUE),
    max_correlation = max(cross_corr, na.rm = TRUE),
    sd_correlation = sd(as.vector(cross_corr), na.rm = TRUE),
    n_comparisons = length(cross_corr),
    high_correlation_count = sum(cross_corr > 0.7, na.rm = TRUE),
    low_correlation_count = sum(cross_corr < 0.3, na.rm = TRUE)
  )
  
  return(stats)
}

# Function to create correlation summary plot
create_correlation_summary_plot <- function(correlation_stats_list, output_file) {
  
  if (length(correlation_stats_list) == 0) {
    cat("No correlation statistics to plot\n")
    return(NULL)
  }
  
  # Prepare data for plotting
  plot_data <- map_dfr(correlation_stats_list, ~ {
    data.frame(
      mean_correlation = .x$mean_correlation,
      median_correlation = .x$median_correlation,
      sd_correlation = .x$sd_correlation,
      high_corr_count = .x$high_correlation_count,
      low_corr_count = .x$low_correlation_count,
      n_comparisons = .x$n_comparisons
    )
  }, .id = "gene_set")
  
  # Create plots
  p1 <- ggplot(plot_data, aes(x = gene_set, y = mean_correlation)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_correlation - sd_correlation, 
                     ymax = mean_correlation + sd_correlation), 
                 width = 0.2) +
    labs(title = "Mean Cross-Dataset Correlation", 
         x = "Gene Set", y = "Mean Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(plot_data, aes(x = gene_set)) +
    geom_col(aes(y = high_corr_count), fill = "green", alpha = 0.7) +
    labs(title = "High Correlation Count (>0.7)", 
         x = "Gene Set", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(plot_data, aes(x = gene_set, y = n_comparisons)) +
    geom_col(fill = "orange", alpha = 0.7) +
    labs(title = "Number of Comparisons", 
         x = "Gene Set", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine plots
  combined_plot <- p1 / p2 / p3
  
  ggsave(output_file, combined_plot, width = 12, height = 10, dpi = 300)
  cat("Correlation summary plot saved:", output_file, "\n")
  
  return(plot_data)
}

# Function to create detailed correlation report
create_correlation_report <- function(correlation_stats_list, plot_data, output_file) {
  
  report_text <- paste0(
    "Cross-Dataset Correlation Analysis Report\n",
    "========================================\n\n",
    "Analysis Date: ", Sys.Date(), "\n",
    "Job ID: ", JOB_ID, "\n",
    "Datasets: ", DATASET_1_NAME, " vs ", DATASET_2_NAME, "\n",
    "Correlation Method: ", CORRELATION_METHOD, "\n\n",
    "Gene Set Analysis Results:\n"
  )
  
  for (gene_set in names(correlation_stats_list)) {
    stats <- correlation_stats_list[[gene_set]]
    report_text <- paste0(report_text,
      "\n", gene_set, ":\n",
      "  - Mean correlation: ", round(stats$mean_correlation, 3), "\n",
      "  - Median correlation: ", round(stats$median_correlation, 3), "\n",
      "  - SD correlation: ", round(stats$sd_correlation, 3), "\n",
      "  - High correlations (>0.7): ", stats$high_correlation_count, "\n",
      "  - Low correlations (<0.3): ", stats$low_correlation_count, "\n",
      "  - Total comparisons: ", stats$n_comparisons, "\n"
    )
  }
  
  # Find best performing gene set
  if (nrow(plot_data) > 0) {
    best_gene_set <- plot_data$gene_set[which.max(plot_data$mean_correlation)]
    report_text <- paste0(report_text,
      "\nBest performing gene set: ", best_gene_set, "\n",
      "Mean correlation: ", round(max(plot_data$mean_correlation), 3), "\n"
    )
  }
  
  report_text <- paste0(report_text,
    "\nRecommendations:\n",
    "1. Focus on gene sets with highest mean correlation\n",
    "2. Investigate biological relevance of highly correlated cell types\n",
    "3. Consider batch effects if correlations are unexpectedly low\n",
    "4. Validate findings with independent datasets\n",
    "5. Perform functional enrichment on highly correlated genes\n\n",
    "Usage Examples:\n",
    "Rscript cross_dataset_correlation.R <input.rds> <dataset1> <dataset2> <job_id>\n",
    "Rscript cross_dataset_correlation.R data.rds CS9 CS10 correlation_v1\n"
  )
  
  writeLines(report_text, output_file)
  cat("Correlation report saved:", output_file, "\n")
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

cat("Starting cross-dataset correlation analysis...\n")
cat("Datasets:", DATASET_1_NAME, "vs", DATASET_2_NAME, "\n")
cat("Job ID:", JOB_ID, "\n")

# Load data
if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded Seurat object with", ncol(sobj), "cells and", nrow(sobj), "features\n")

# Validate metadata columns
if (!CONDITION_COL %in% colnames(sobj@meta.data)) {
  stop("Condition column not found: ", CONDITION_COL)
}

if (!CELL_TYPE_COL %in% colnames(sobj@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    CELL_TYPE_COL <- "preliminary_celltype"
  } else {
    stop("Cell type column not found")
  }
}

# Check datasets
available_conditions <- unique(sobj@meta.data[[CONDITION_COL]])
cat("Available conditions:", paste(available_conditions, collapse = ", "), "\n")

if (!DATASET_1_NAME %in% available_conditions) {
  stop("Dataset 1 not found: ", DATASET_1_NAME)
}
if (!DATASET_2_NAME %in% available_conditions) {
  stop("Dataset 2 not found: ", DATASET_2_NAME)
}

# ============================================================================
# Prepare Expression Matrices
# ============================================================================

cat("Preparing expression matrices...\n")

# Subset data by condition
sobj_dataset1 <- subset(sobj, cells = which(sobj@meta.data[[CONDITION_COL]] == DATASET_1_NAME))
sobj_dataset2 <- subset(sobj, cells = which(sobj@meta.data[[CONDITION_COL]] == DATASET_2_NAME))

cat("Dataset 1 (", DATASET_1_NAME, "):", ncol(sobj_dataset1), "cells\n")
cat("Dataset 2 (", DATASET_2_NAME, "):", ncol(sobj_dataset2), "cells\n")

# Prepare expression matrices
mat1 <- prepare_expression_matrix(sobj_dataset1, DATASET_1_NAME)
mat2 <- prepare_expression_matrix(sobj_dataset2, DATASET_2_NAME)

cat("Expression matrix 1:", nrow(mat1), "genes x", ncol(mat1), "cell types\n")
cat("Expression matrix 2:", nrow(mat2), "genes x", ncol(mat2), "cell types\n")

# Find common genes
all_common_genes <- intersect(rownames(mat1), rownames(mat2))
cat("Common genes between datasets:", length(all_common_genes), "\n")

if (length(all_common_genes) == 0) {
  stop("No common genes found between datasets")
}

# ============================================================================
# Load DEG Data for Gene Selection
# ============================================================================

cat("Loading DEG data for gene selection...\n")

# Try to load DEG data
deg_file_1 <- file.path(DATA_DIR, paste0("DEG_", DATASET_1_NAME, "_", JOB_ID, ".rds"))
deg_file_2 <- file.path(DATA_DIR, paste0("DEG_", DATASET_2_NAME, "_", JOB_ID, ".rds"))

deg_data_1 <- load_deg_data(deg_file_1)
deg_data_2 <- load_deg_data(deg_file_2)

# If DEG data not available, create mock data or use all genes
if (is.null(deg_data_1) && is.null(deg_data_2)) {
  cat("No DEG data found, using all common genes\n")
  gene_sets <- list("all_genes" = all_common_genes)
} else {
  # Create gene sets based on DEG data
  gene_sets <- list()
  
  for (top_n in TOP_GENES_OPTIONS) {
    genes_1 <- if (!is.null(deg_data_1)) select_top_genes(deg_data_1, top_n) else character(0)
    genes_2 <- if (!is.null(deg_data_2)) select_top_genes(deg_data_2, top_n) else character(0)
    
    # Find common genes among DEGs
    if (length(genes_1) > 0 && length(genes_2) > 0) {
      common_deg_genes <- intersect(intersect(genes_1, genes_2), all_common_genes)
    } else if (length(genes_1) > 0) {
      common_deg_genes <- intersect(genes_1, all_common_genes)
    } else if (length(genes_2) > 0) {
      common_deg_genes <- intersect(genes_2, all_common_genes)
    } else {
      common_deg_genes <- character(0)
    }
    
    if (length(common_deg_genes) > 0) {
      gene_sets[[paste0("top_", top_n, "_genes")]] <- common_deg_genes
    }
  }
  
  # Always include all common genes as a reference
  gene_sets[["all_genes"]] <- all_common_genes
}

cat("Gene sets created:", length(gene_sets), "\n")
for (set_name in names(gene_sets)) {
  cat("  -", set_name, ":", length(gene_sets[[set_name]]), "genes\n")
}

# ============================================================================
# Perform Correlation Analysis
# ============================================================================

cat("Performing correlation analysis...\n")

correlation_matrices <- list()
correlation_stats_list <- list()

for (set_name in names(gene_sets)) {
  cat("Analyzing gene set:", set_name, "\n")
  
  genes <- gene_sets[[set_name]]
  
  # Calculate correlation matrix
  corr_matrix <- calculate_correlation_matrix(mat1, mat2, genes)
  
  if (!is.null(corr_matrix)) {
    correlation_matrices[[set_name]] <- corr_matrix
    
    # Calculate statistics
    stats <- calculate_correlation_stats(corr_matrix, DATASET_1_NAME, DATASET_2_NAME)
    correlation_stats_list[[set_name]] <- stats
    
    # Create heatmap
    output_file <- file.path(PLOTS_DIR, 
                           paste0("correlation_heatmap_", set_name, "_", JOB_ID, ".pdf"))
    title <- paste("Correlation Heatmap -", set_name, 
                  paste0("(", length(genes), " genes)"))
    
    create_correlation_heatmap(corr_matrix, title, output_file)
  }
}

# ============================================================================
# Create Summary Visualizations and Reports
# ============================================================================

cat("Creating summary visualizations...\n")

# Create correlation summary plot
if (length(correlation_stats_list) > 0) {
  summary_plot_file <- file.path(PLOTS_DIR, paste0("correlation_summary_", JOB_ID, ".pdf"))
  plot_data <- create_correlation_summary_plot(correlation_stats_list, summary_plot_file)
  
  # Save correlation statistics
  stats_file <- file.path(OUTPUT_DIR, paste0("correlation_statistics_", JOB_ID, ".xlsx"))
  
  # Prepare data for Excel export
  excel_data <- list()
  
  if (!is.null(plot_data)) {
    excel_data[["Summary"]] <- plot_data
  }
  
  for (set_name in names(correlation_stats_list)) {
    stats_df <- data.frame(
      metric = names(correlation_stats_list[[set_name]]),
      value = unlist(correlation_stats_list[[set_name]])
    )
    excel_data[[paste0("Stats_", set_name)]] <- stats_df
  }
  
  if (length(excel_data) > 0) {
    write_xlsx(excel_data, stats_file)
    cat("Correlation statistics saved:", stats_file, "\n")
  }
}

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating comprehensive report...\n")

report_file <- file.path(OUTPUT_DIR, paste0("correlation_analysis_report_", JOB_ID, ".txt"))
create_correlation_report(correlation_stats_list, plot_data, report_file)

# Save session info
session_file <- file.path(OUTPUT_DIR, paste0("correlation_sessionInfo_", JOB_ID, ".txt"))
writeLines(capture.output(sessionInfo()), session_file)

# ============================================================================
# Final Summary
# ============================================================================

cat("\n", "="*60, "\n")
cat("Cross-dataset correlation analysis completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\nFinal Summary:\n")
cat("- Datasets compared:", DATASET_1_NAME, "vs", DATASET_2_NAME, "\n")
cat("- Common genes:", length(all_common_genes), "\n")
cat("- Gene sets analyzed:", length(gene_sets), "\n")
cat("- Correlation method:", CORRELATION_METHOD, "\n")

if (length(correlation_stats_list) > 0) {
  best_set <- names(correlation_stats_list)[which.max(sapply(correlation_stats_list, function(x) x$mean_correlation))]
  best_corr <- correlation_stats_list[[best_set]]$mean_correlation
  cat("- Best gene set:", best_set, "(mean correlation:", round(best_corr, 3), ")\n")
}

cat("\nNext steps:\n")
cat("1. Review correlation heatmaps for biological insights\n")
cat("2. Investigate highly correlated cell type pairs\n")
cat("3. Validate findings with independent data\n")
cat("4. Consider functional enrichment analysis\n")

cat("\nAnalysis completed successfully!\n")