#!/usr/bin/env Rscript

# ============================================================================
# Correlation Analysis for Integrated Datasets
# ============================================================================
# 
# This script performs comprehensive correlation analysis between different
# datasets, developmental stages, or cell types. It generates correlation
# heatmaps and quantitative metrics for comparative analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(corrplot)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
})

# Parameters
CORRELATION_METHOD <- "spearman"    # Correlation method: "pearson" or "spearman"
TOP_GENES_OPTIONS <- c(25, 50, 100, 200, 500)  # Number of top DEGs to use
CLUSTER_HEATMAP <- TRUE            # Whether to cluster rows/columns in heatmap
MIN_CELLS_PER_TYPE <- 10           # Minimum cells per cell type for inclusion

# Input/Output paths
INPUT_SEURAT1 <- "../data/dataset1_with_annotations.rds"  # First dataset
INPUT_SEURAT2 <- "../data/dataset2_with_annotations.rds"  # Second dataset (optional)
INPUT_DEG1 <- "../data/deg_markers_dataset1.rds"         # DEG markers for dataset 1
INPUT_DEG2 <- "../data/deg_markers_dataset2.rds"         # DEG markers for dataset 2 (optional)
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Helper Functions
# ============================================================================

# Function to prepare expression matrix from Seurat object
prepare_expression_matrix <- function(sobj, group_by = "celltype") {
  cat("Preparing expression matrix for", ncol(sobj), "cells...\\n")
  
  # Filter cell types with sufficient cells
  cell_counts <- table(sobj@meta.data[[group_by]])
  valid_types <- names(cell_counts)[cell_counts >= MIN_CELLS_PER_TYPE]
  
  if (length(valid_types) == 0) {
    stop("No cell types have sufficient cells (>= ", MIN_CELLS_PER_TYPE, ")")
  }
  
  # Subset to valid cell types
  sobj_filtered <- subset(sobj, subset = eval(parse(text = paste0(group_by, " %in% c('", 
                                                                  paste(valid_types, collapse = "', '"), "')"))))
  
  cat("Using", length(valid_types), "cell types with >= ", MIN_CELLS_PER_TYPE, "cells\\n")
  
  # Aggregate expression by cell type
  Idents(sobj_filtered) <- group_by
  avgexp <- AverageExpression(sobj_filtered, return.seurat = FALSE, assays = "RNA")
  
  return(avgexp$RNA)
}

# Function to select top DEG genes
select_top_deg_genes <- function(deg_data, top_n = 100, method = "avg_log2FC") {
  if (is.null(deg_data) || nrow(deg_data) == 0) {
    return(character(0))
  }
  
  top_genes <- deg_data %>%
    group_by(cluster) %>%
    arrange(desc(!!sym(method))) %>%
    slice_head(n = top_n) %>%
    pull(gene) %>%
    unique()
  
  cat("Selected", length(top_genes), "top DEG genes\\n")
  return(top_genes)
}

# Function to create correlation heatmap
create_correlation_heatmap <- function(cor_matrix, title = "Correlation Heatmap", 
                                     filename = NULL, cluster = TRUE) {
  
  # Define color palette
  colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  breaks <- seq(-1, 1, length.out = 101)
  
  # Create heatmap
  if (!is.null(filename)) {
    pdf(filename, width = 12, height = 10)
  }
  
  p <- pheatmap(
    cor_matrix,
    cluster_rows = cluster,
    cluster_cols = cluster,
    color = colors,
    breaks = breaks,
    display_numbers = TRUE,
    number_color = "black",
    number_format = "%.2f",
    fontsize = 8,
    fontsize_number = 6,
    main = title,
    border_color = "white",
    cellwidth = 15,
    cellheight = 15
  )
  
  if (!is.null(filename)) {
    dev.off()
    cat("Heatmap saved to:", filename, "\\n")
  }
  
  return(p)
}

# Function to calculate correlation statistics
calculate_correlation_stats <- function(cor_matrix) {
  # Remove diagonal (self-correlations)
  diag(cor_matrix) <- NA
  
  stats <- list(
    mean_correlation = mean(cor_matrix, na.rm = TRUE),
    median_correlation = median(cor_matrix, na.rm = TRUE),
    min_correlation = min(cor_matrix, na.rm = TRUE),
    max_correlation = max(cor_matrix, na.rm = TRUE),
    sd_correlation = sd(cor_matrix, na.rm = TRUE),
    high_correlations = sum(cor_matrix > 0.7, na.rm = TRUE),
    moderate_correlations = sum(cor_matrix > 0.5 & cor_matrix <= 0.7, na.rm = TRUE),
    low_correlations = sum(cor_matrix <= 0.5, na.rm = TRUE)
  )
  
  return(stats)
}

# ============================================================================
# Load Data
# ============================================================================

cat("Loading datasets...\\n")

# Load first dataset
if (!file.exists(INPUT_SEURAT1)) {
  stop("First dataset not found: ", INPUT_SEURAT1)
}
sobj1 <- readRDS(INPUT_SEURAT1)
cat("Loaded dataset 1 with", ncol(sobj1), "cells\\n")

# Load second dataset (if provided)
sobj2 <- NULL
if (file.exists(INPUT_SEURAT2)) {
  sobj2 <- readRDS(INPUT_SEURAT2)
  cat("Loaded dataset 2 with", ncol(sobj2), "cells\\n")
} else {
  cat("Second dataset not found, performing single-dataset correlation analysis\\n")
}

# Load DEG data
deg1 <- NULL
if (file.exists(INPUT_DEG1)) {
  deg1 <- readRDS(INPUT_DEG1)
  cat("Loaded DEG data for dataset 1 with", nrow(deg1), "markers\\n")
}

deg2 <- NULL
if (!is.null(sobj2) && file.exists(INPUT_DEG2)) {
  deg2 <- readRDS(INPUT_DEG2)
  cat("Loaded DEG data for dataset 2 with", nrow(deg2), "markers\\n")
}

# ============================================================================
# Prepare Expression Data
# ============================================================================

cat("\\nPreparing expression matrices...\\n")

# Determine grouping variable
group_by_var <- "celltype"
if (!"celltype" %in% colnames(sobj1@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj1@meta.data)) {
    group_by_var <- "preliminary_celltype"
  } else if ("seurat_clusters" %in% colnames(sobj1@meta.data)) {
    group_by_var <- "seurat_clusters"
  } else {
    stop("No suitable grouping variable found (celltype, preliminary_celltype, or seurat_clusters)")
  }
}

cat("Using grouping variable:", group_by_var, "\\n")

# Prepare expression matrices
expr_matrix1 <- prepare_expression_matrix(sobj1, group_by = group_by_var)

expr_matrix2 <- NULL
if (!is.null(sobj2)) {
  expr_matrix2 <- prepare_expression_matrix(sobj2, group_by = group_by_var)
}

# ============================================================================
# Correlation Analysis
# ============================================================================

cat("\\nPerforming correlation analysis...\\n")

correlation_results <- list()

# Single dataset analysis (intra-dataset correlations)
if (is.null(sobj2)) {
  cat("Performing intra-dataset correlation analysis...\\n")
  
  # All genes correlation
  cor_all_genes <- cor(expr_matrix1, method = CORRELATION_METHOD)
  correlation_results[["all_genes"]] <- cor_all_genes
  
  # Create heatmap for all genes
  create_correlation_heatmap(
    cor_all_genes,
    title = paste("Intra-dataset Correlation (All Genes,", CORRELATION_METHOD, ")"),
    filename = file.path(PLOTS_DIR, "correlation_heatmap_all_genes.pdf"),
    cluster = CLUSTER_HEATMAP
  )
  
  # DEG-based correlations
  if (!is.null(deg1)) {
    for (top_n in TOP_GENES_OPTIONS) {
      cat("Analyzing with top", top_n, "DEGs...\\n")
      
      top_genes <- select_top_deg_genes(deg1, top_n = top_n)
      
      if (length(top_genes) > 0) {
        # Filter expression matrix to top genes
        common_genes <- intersect(top_genes, rownames(expr_matrix1))
        
        if (length(common_genes) >= 10) {
          expr_filtered <- expr_matrix1[common_genes, ]
          cor_deg <- cor(expr_filtered, method = CORRELATION_METHOD)
          correlation_results[[paste0("top_", top_n, "_degs")]] <- cor_deg
          
          # Create heatmap
          create_correlation_heatmap(
            cor_deg,
            title = paste("Intra-dataset Correlation (Top", top_n, "DEGs,", CORRELATION_METHOD, ")"),
            filename = file.path(PLOTS_DIR, paste0("correlation_heatmap_top_", top_n, "_degs.pdf")),
            cluster = CLUSTER_HEATMAP
          )
        } else {
          cat("Warning: Only", length(common_genes), "common genes found for top", top_n, "DEGs\\n")
        }
      }
    }
  }
  
} else {
  # Cross-dataset analysis (inter-dataset correlations)
  cat("Performing inter-dataset correlation analysis...\\n")
  
  # Find common cell types
  types1 <- colnames(expr_matrix1)
  types2 <- colnames(expr_matrix2)
  
  cat("Dataset 1 cell types:", length(types1), "\\n")
  cat("Dataset 2 cell types:", length(types2), "\\n")
  
  # All genes correlation
  common_genes_all <- intersect(rownames(expr_matrix1), rownames(expr_matrix2))
  cat("Common genes between datasets:", length(common_genes_all), "\\n")
  
  if (length(common_genes_all) >= 100) {
    expr1_filtered <- expr_matrix1[common_genes_all, ]
    expr2_filtered <- expr_matrix2[common_genes_all, ]
    
    # Combine matrices
    combined_matrix <- cbind(expr1_filtered, expr2_filtered)
    colnames(combined_matrix) <- c(paste0("Dataset1_", colnames(expr1_filtered)),
                                  paste0("Dataset2_", colnames(expr2_filtered)))
    
    cor_cross <- cor(combined_matrix, method = CORRELATION_METHOD)
    correlation_results[["cross_dataset_all_genes"]] <- cor_cross
    
    # Create heatmap
    create_correlation_heatmap(
      cor_cross,
      title = paste("Cross-dataset Correlation (All Genes,", CORRELATION_METHOD, ")"),
      filename = file.path(PLOTS_DIR, "correlation_heatmap_cross_dataset_all_genes.pdf"),
      cluster = CLUSTER_HEATMAP
    )
    
    # Extract cross-correlations only (dataset1 vs dataset2)
    cross_only <- cor_cross[1:ncol(expr1_filtered), 
                           (ncol(expr1_filtered)+1):ncol(combined_matrix)]
    correlation_results[["cross_correlations_only"]] <- cross_only
    
    create_correlation_heatmap(
      cross_only,
      title = paste("Dataset1 vs Dataset2 Correlation (All Genes,", CORRELATION_METHOD, ")"),
      filename = file.path(PLOTS_DIR, "correlation_heatmap_cross_only_all_genes.pdf"),
      cluster = CLUSTER_HEATMAP
    )
  }
  
  # DEG-based cross-dataset correlations
  if (!is.null(deg1) && !is.null(deg2)) {
    for (top_n in TOP_GENES_OPTIONS) {
      cat("Cross-dataset analysis with top", top_n, "DEGs...\\n")
      
      top_genes1 <- select_top_deg_genes(deg1, top_n = top_n)
      top_genes2 <- select_top_deg_genes(deg2, top_n = top_n)
      
      # Use union of top genes from both datasets
      union_genes <- union(top_genes1, top_genes2)
      common_deg_genes <- intersect(union_genes, common_genes_all)
      
      if (length(common_deg_genes) >= 10) {
        expr1_deg <- expr_matrix1[common_deg_genes, ]
        expr2_deg <- expr_matrix2[common_deg_genes, ]
        
        combined_deg_matrix <- cbind(expr1_deg, expr2_deg)
        colnames(combined_deg_matrix) <- c(paste0("Dataset1_", colnames(expr1_deg)),
                                          paste0("Dataset2_", colnames(expr2_deg)))
        
        cor_cross_deg <- cor(combined_deg_matrix, method = CORRELATION_METHOD)
        correlation_results[[paste0("cross_dataset_top_", top_n, "_degs")]] <- cor_cross_deg
        
        # Create heatmap
        create_correlation_heatmap(
          cor_cross_deg,
          title = paste("Cross-dataset Correlation (Top", top_n, "DEGs,", CORRELATION_METHOD, ")"),
          filename = file.path(PLOTS_DIR, paste0("correlation_heatmap_cross_dataset_top_", top_n, "_degs.pdf")),
          cluster = CLUSTER_HEATMAP
        )
        
        # Cross-correlations only
        cross_deg_only <- cor_cross_deg[1:ncol(expr1_deg), 
                                       (ncol(expr1_deg)+1):ncol(combined_deg_matrix)]
        correlation_results[[paste0("cross_correlations_top_", top_n, "_degs")]] <- cross_deg_only
        
        create_correlation_heatmap(
          cross_deg_only,
          title = paste("Dataset1 vs Dataset2 (Top", top_n, "DEGs,", CORRELATION_METHOD, ")"),
          filename = file.path(PLOTS_DIR, paste0("correlation_heatmap_cross_only_top_", top_n, "_degs.pdf")),
          cluster = CLUSTER_HEATMAP
        )
      } else {
        cat("Warning: Only", length(common_deg_genes), "common DEG genes found for top", top_n, "\\n")
      }
    }
  }
}

# ============================================================================
# Statistical Analysis
# ============================================================================

cat("\\nCalculating correlation statistics...\\n")

# Calculate statistics for all correlation matrices
correlation_stats <- list()
for (analysis_name in names(correlation_results)) {
  cor_matrix <- correlation_results[[analysis_name]]
  stats <- calculate_correlation_stats(cor_matrix)
  correlation_stats[[analysis_name]] <- stats
  
  cat("\\n", analysis_name, "statistics:\\n")
  cat("  Mean correlation:", round(stats$mean_correlation, 3), "\\n")
  cat("  High correlations (>0.7):", stats$high_correlations, "\\n")
  cat("  Moderate correlations (0.5-0.7):", stats$moderate_correlations, "\\n")
  cat("  Low correlations (<=0.5):", stats$low_correlations, "\\n")
}

# Convert statistics to data frame
stats_df <- map_dfr(correlation_stats, ~ as.data.frame(.x), .id = "analysis")

# ============================================================================
# Generate Summary Report
# ============================================================================

cat("\\nGenerating correlation analysis report...\\n")

# Create summary statistics
n_analyses <- length(correlation_results)
analysis_type <- ifelse(is.null(sobj2), "Intra-dataset", "Cross-dataset")

# Generate report
report_text <- paste0(
  "Correlation Analysis Report\\n",
  "=========================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Analysis Type: ", analysis_type, "\\n",
  "Correlation Method: ", toupper(CORRELATION_METHOD), "\\n\\n",
  "Dataset Information:\\n",
  "  - Dataset 1: ", ncol(sobj1), " cells, ", nrow(expr_matrix1), " genes\\n"
)

if (!is.null(sobj2)) {
  report_text <- paste0(report_text,
    "  - Dataset 2: ", ncol(sobj2), " cells, ", nrow(expr_matrix2), " genes\\n"
  )
}

report_text <- paste0(report_text,
  "  - Grouping variable: ", group_by_var, "\\n",
  "  - Minimum cells per type: ", MIN_CELLS_PER_TYPE, "\\n\\n",
  "Analysis Parameters:\\n",
  "  - Number of analyses performed: ", n_analyses, "\\n",
  "  - Top DEG options tested: ", paste(TOP_GENES_OPTIONS, collapse = ", "), "\\n",
  "  - Heatmap clustering: ", CLUSTER_HEATMAP, "\\n\\n",
  "Results Summary:\\n"
)

# Add statistics for each analysis
for (analysis_name in names(correlation_stats)) {
  stats <- correlation_stats[[analysis_name]]
  report_text <- paste0(report_text,
    "\\n", analysis_name, ":\\n",
    "  - Mean correlation: ", round(stats$mean_correlation, 3), "\\n",
    "  - Range: [", round(stats$min_correlation, 3), ", ", 
         round(stats$max_correlation, 3), "]\\n",
    "  - High correlations (>0.7): ", stats$high_correlations, "\\n"
  )
}

report_text <- paste0(report_text,
  "\\nOutput Files:\\n",
  "  - Correlation matrices: correlation_results.rds\\n",
  "  - Statistics summary: correlation_statistics.csv\\n",
  "  - Individual heatmaps: correlation_heatmap_*.pdf\\n\\n",
  "Recommendations:\\n",
  "1. Review correlation patterns for biological consistency\\n",
  "2. Investigate highly correlated cell types for potential merging\\n",
  "3. Validate low correlations with additional markers\\n",
  "4. Consider different gene selection strategies if needed\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "correlation_analysis_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("\\nSaving results...\\n")

# Save correlation matrices
saveRDS(correlation_results, file.path(OUTPUT_DIR, "correlation_results.rds"))

# Save statistics
write.csv(stats_df, file.path(OUTPUT_DIR, "correlation_statistics.csv"), 
          row.names = FALSE)

# Save processed expression matrices
saveRDS(expr_matrix1, file.path(OUTPUT_DIR, "expression_matrix_dataset1.rds"))
if (!is.null(expr_matrix2)) {
  saveRDS(expr_matrix2, file.path(OUTPUT_DIR, "expression_matrix_dataset2.rds"))
}

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_correlation_analysis.txt"))

cat("\\nCorrelation analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nSummary:\\n")
cat("- Performed", n_analyses, "correlation analyses\\n")
cat("- Generated", length(list.files(PLOTS_DIR, pattern = "correlation_heatmap.*\\\\.pdf")), "heatmap visualizations\\n")
cat("- Analysis type:", analysis_type, "\\n")
cat("\\nNext steps:\\n")
cat("1. Review heatmaps for biological patterns\\n")
cat("2. Investigate unexpected correlations\\n") 
cat("3. Use results for downstream comparative analysis\\n")