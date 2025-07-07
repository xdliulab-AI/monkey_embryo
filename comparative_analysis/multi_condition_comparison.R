#!/usr/bin/env Rscript

# ============================================================================
# Multi-Condition Comparative Analysis
# ============================================================================
# 
# This script performs systematic comparative analysis across multiple 
# conditions, datasets, or cell lineages. It includes batch processing
# of multiple Seurat objects, integration analysis, and comprehensive
# comparison workflows.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(openxlsx)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(ComplexHeatmap)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(future.apply)
  library(gridExtra)
  library(corrplot)
})

# Parameters
JOB_ID <- "multi_comparison"
CONDITIONS <- c("CS9", "CS10")         # Conditions to compare
CELL_TYPE_COL <- "celltype"           # Column name for cell types
CONDITION_COL <- "CS"                 # Column name for conditions
MIN_PCT <- 0.25                       # Minimum percentage for DEG analysis
LOGFC_THRESHOLD <- 0.25               # Log fold change threshold
P_VAL_CUTOFF <- 0.05                  # P-value cutoff
N_WORKERS <- 4                        # Number of parallel workers
BATCH_PROCESS <- TRUE                 # Process multiple files

# File processing parameters
FILE_PATTERNS <- c(
  "*Paraxial*Mesoderm*.rds",
  "*Brain*spinal*cord*.rds", 
  "*Lateral*Plate*Mesoderm*.rds"
)

# Input/Output paths
INPUT_DIR <- "../data"
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
  INPUT_DIR <- args[1]
}
if (length(args) > 1) {
  JOB_ID <- args[2]
}
if (length(args) > 2) {
  CONDITIONS <- strsplit(args[3], ",")[[1]]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to find input files
find_input_files <- function(input_dir, patterns = FILE_PATTERNS) {
  all_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(patterns) > 0) {
    matched_files <- c()
    for (pattern in patterns) {
      pattern_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
      matched_files <- c(matched_files, pattern_files)
    }
    matched_files <- unique(matched_files)
  } else {
    matched_files <- all_files
  }
  
  return(matched_files)
}

# Function to extract lineage name from filename
extract_lineage_name <- function(filename) {
  base_name <- basename(filename)
  base_name <- gsub("\\.rds$", "", base_name)
  base_name <- gsub("^[0-9_]*", "", base_name)  # Remove leading numbers
  base_name <- gsub("sobj_", "", base_name)     # Remove sobj prefix
  base_name <- gsub("_res[0-9.]+.*", "", base_name)  # Remove resolution suffix
  return(base_name)
}

# Function to validate Seurat object
validate_seurat_object <- function(sobj, filename) {
  errors <- c()
  
  if (is.null(sobj)) {
    errors <- c(errors, "Object is NULL")
    return(errors)
  }
  
  if (!inherits(sobj, "Seurat")) {
    errors <- c(errors, "Not a Seurat object")
    return(errors)
  }
  
  if (ncol(sobj) == 0) {
    errors <- c(errors, "No cells in object")
  }
  
  if (nrow(sobj) == 0) {
    errors <- c(errors, "No genes in object")
  }
  
  if (!CONDITION_COL %in% colnames(sobj@meta.data)) {
    errors <- c(errors, paste("Missing condition column:", CONDITION_COL))
  }
  
  if (!CELL_TYPE_COL %in% colnames(sobj@meta.data)) {
    if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
      # This is acceptable, we'll use it instead
    } else {
      errors <- c(errors, paste("Missing cell type column:", CELL_TYPE_COL))
    }
  }
  
  if (length(errors) > 0) {
    cat("Validation errors for", basename(filename), ":\n")
    cat(paste("  -", errors, collapse = "\n"), "\n")
  }
  
  return(errors)
}

# Function to perform DEG analysis for one condition pair
perform_deg_analysis <- function(sobj, cond1, cond2, lineage_name) {
  
  # Set cell type as identity
  cell_type_col <- if (CELL_TYPE_COL %in% colnames(sobj@meta.data)) {
    CELL_TYPE_COL
  } else if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    "preliminary_celltype"
  } else {
    "seurat_clusters"
  }
  
  Idents(sobj) <- cell_type_col
  
  # Get available conditions
  available_conditions <- unique(sobj@meta.data[[CONDITION_COL]])
  
  if (!all(c(cond1, cond2) %in% available_conditions)) {
    missing_conds <- setdiff(c(cond1, cond2), available_conditions)
    cat("Missing conditions in", lineage_name, ":", paste(missing_conds, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Get clusters and sizes
  clusters <- Idents(sobj) %>% unique() %>% sort()
  cluster_sizes <- table(Idents(sobj))
  
  # Function for single cluster analysis
  analyze_cluster <- function(cluster_id) {
    if (cluster_sizes[[cluster_id]] < 3) {
      return(NULL)
    }
    
    sobj_cluster <- subset(sobj, idents = cluster_id)
    
    # Check if both conditions present
    cluster_conditions <- unique(sobj_cluster@meta.data[[CONDITION_COL]])
    if (!all(c(cond1, cond2) %in% cluster_conditions)) {
      return(NULL)
    }
    
    # Check cell counts per condition
    condition_counts <- table(sobj_cluster@meta.data[[CONDITION_COL]])
    if (any(condition_counts < 3)) {
      return(NULL)
    }
    
    # Perform DEG analysis
    tryCatch({
      markers <- FindMarkers(
        sobj_cluster,
        ident.1 = cond1,
        ident.2 = cond2,
        group.by = CONDITION_COL,
        only.pos = FALSE,
        min.pct = MIN_PCT,
        logfc.threshold = LOGFC_THRESHOLD,
        test.use = "wilcox"
      )
      
      if (nrow(markers) > 0) {
        markers <- markers %>%
          as.data.frame() %>%
          mutate(
            gene = rownames(.),
            cluster = cluster_id,
            lineage = lineage_name,
            comparison = paste(cond1, "vs", cond2),
            significant = ifelse(p_val_adj < P_VAL_CUTOFF, "Significant", "Not significant"),
            regulation = ifelse(avg_log2FC > 0, paste("Up in", cond1), paste("Up in", cond2)),
            abs_log2FC = abs(avg_log2FC)
          ) %>%
          select(gene, cluster, lineage, comparison, everything()) %>%
          arrange(desc(abs_log2FC))
        
        return(markers)
      } else {
        return(NULL)
      }
    }, error = function(e) {
      cat("Error analyzing cluster", cluster_id, "in", lineage_name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Analyze all clusters
  cluster_results <- map(clusters, analyze_cluster)
  names(cluster_results) <- clusters
  cluster_results <- discard(cluster_results, is.null)
  
  if (length(cluster_results) > 0) {
    combined_results <- bind_rows(cluster_results)
    cat("Found", nrow(combined_results), "DEGs in", lineage_name, "\n")
    return(combined_results)
  } else {
    cat("No DEGs found in", lineage_name, "\n")
    return(NULL)
  }
}

# Function to create comparison summary
create_comparison_summary <- function(all_results) {
  
  if (is.null(all_results) || nrow(all_results) == 0) {
    return(NULL)
  }
  
  # Summary by lineage
  lineage_summary <- all_results %>%
    group_by(lineage, comparison) %>%
    summarise(
      total_genes = n(),
      significant_genes = sum(significant == "Significant"),
      up_regulated = sum(str_detect(regulation, paste("Up in", CONDITIONS[1]))),
      down_regulated = sum(str_detect(regulation, paste("Up in", CONDITIONS[2]))),
      mean_log2FC = mean(abs_log2FC),
      clusters_analyzed = n_distinct(cluster),
      .groups = "drop"
    ) %>%
    arrange(desc(total_genes))
  
  # Summary by cluster across lineages
  cluster_summary <- all_results %>%
    group_by(cluster, lineage) %>%
    summarise(
      total_genes = n(),
      significant_genes = sum(significant == "Significant"),
      mean_log2FC = mean(abs_log2FC),
      .groups = "drop"
    ) %>%
    arrange(desc(total_genes))
  
  # Overall summary
  overall_summary <- all_results %>%
    summarise(
      total_lineages = n_distinct(lineage),
      total_clusters = n_distinct(paste(lineage, cluster)),
      total_genes = n(),
      significant_genes = sum(significant == "Significant"),
      mean_log2FC = mean(abs_log2FC),
      .groups = "drop"
    )
  
  return(list(
    lineage = lineage_summary,
    cluster = cluster_summary,
    overall = overall_summary
  ))
}

# Function to create cross-lineage comparison plots
create_cross_lineage_plots <- function(all_results, output_prefix) {
  
  if (is.null(all_results) || nrow(all_results) == 0) {
    return(NULL)
  }
  
  plots_created <- c()
  
  # 1. DEG count by lineage
  lineage_counts <- all_results %>%
    group_by(lineage) %>%
    summarise(
      total = n(),
      significant = sum(significant == "Significant"),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(total, significant), names_to = "type", values_to = "count")
  
  p1 <- ggplot(lineage_counts, aes(x = reorder(lineage, count), y = count, fill = type)) +
    geom_col(position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("total" = "lightblue", "significant" = "darkblue")) +
    labs(title = "DEG Count by Lineage", x = "Lineage", y = "Number of DEGs") +
    theme_minimal()
  
  plot_file1 <- file.path(PLOTS_DIR, paste0(output_prefix, "_lineage_deg_counts.pdf"))
  ggsave(plot_file1, p1, width = 12, height = 8, dpi = 300)
  plots_created <- c(plots_created, "lineage_deg_counts")
  
  # 2. Log2FC distribution by lineage
  p2 <- ggplot(all_results, aes(x = lineage, y = avg_log2FC, fill = lineage)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = "Log2FC Distribution by Lineage", 
         x = "Lineage", y = "Average Log2 Fold Change") +
    theme_minimal() +
    theme(legend.position = "none")
  
  plot_file2 <- file.path(PLOTS_DIR, paste0(output_prefix, "_lineage_log2fc_dist.pdf"))
  ggsave(plot_file2, p2, width = 12, height = 8, dpi = 300)
  plots_created <- c(plots_created, "lineage_log2fc_dist")
  
  # 3. Regulation direction by lineage
  regulation_summary <- all_results %>%
    group_by(lineage, regulation) %>%
    summarise(count = n(), .groups = "drop")
  
  p3 <- ggplot(regulation_summary, aes(x = lineage, y = count, fill = regulation)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_fill_manual(values = c("red", "blue")) +
    labs(title = "Regulation Direction by Lineage", 
         x = "Lineage", y = "Number of DEGs") +
    theme_minimal()
  
  plot_file3 <- file.path(PLOTS_DIR, paste0(output_prefix, "_lineage_regulation.pdf"))
  ggsave(plot_file3, p3, width = 12, height = 8, dpi = 300)
  plots_created <- c(plots_created, "lineage_regulation")
  
  return(plots_created)
}

# Function to save comprehensive results
save_comprehensive_results <- function(all_results, summary_data, output_prefix) {
  
  files_saved <- c()
  
  # Save main results
  if (!is.null(all_results) && nrow(all_results) > 0) {
    csv_file <- file.path(DATA_DIR, paste0(output_prefix, "_all_results.csv"))
    rds_file <- file.path(DATA_DIR, paste0(output_prefix, "_all_results.rds"))
    
    write.csv(all_results, csv_file, row.names = FALSE)
    saveRDS(all_results, rds_file)
    
    files_saved <- c(files_saved, csv_file, rds_file)
  }
  
  # Save Excel workbook with multiple sheets
  if (!is.null(summary_data)) {
    xlsx_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_comprehensive_results.xlsx"))
    
    excel_data <- list()
    
    if (!is.null(all_results) && nrow(all_results) > 0) {
      excel_data[["All_Results"]] <- all_results
    }
    
    if (!is.null(summary_data$overall)) {
      excel_data[["Overall_Summary"]] <- summary_data$overall
    }
    
    if (!is.null(summary_data$lineage)) {
      excel_data[["Lineage_Summary"]] <- summary_data$lineage
    }
    
    if (!is.null(summary_data$cluster)) {
      excel_data[["Cluster_Summary"]] <- summary_data$cluster
    }
    
    # Create workbook
    wb <- createWorkbook()
    
    for (sheet_name in names(excel_data)) {
      addWorksheet(wb, sheet_name, zoom = 200)
      writeData(wb, sheet_name, excel_data[[sheet_name]])
    }
    
    saveWorkbook(wb, xlsx_file, overwrite = TRUE)
    files_saved <- c(files_saved, xlsx_file)
  }
  
  return(files_saved)
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

cat("Starting multi-condition comparative analysis...\n")
cat("Job ID:", JOB_ID, "\n")
cat("Conditions:", paste(CONDITIONS, collapse = " vs "), "\n")

# Set up parallel processing
plan(multisession, workers = N_WORKERS)
options(future.globals.maxSize = 3 * 1024^3)

# Find input files
if (BATCH_PROCESS) {
  input_files <- find_input_files(INPUT_DIR, FILE_PATTERNS)
} else {
  input_files <- list.files(INPUT_DIR, pattern = "\\.rds$", full.names = TRUE)
}

cat("Found", length(input_files), "input files to process\n")

if (length(input_files) == 0) {
  stop("No input files found in directory: ", INPUT_DIR)
}

# Process each file
all_lineage_results <- list()

for (file_path in input_files) {
  lineage_name <- extract_lineage_name(file_path)
  cat("\nProcessing:", lineage_name, "\n")
  cat("File:", basename(file_path), "\n")
  
  # Load Seurat object
  tryCatch({
    sobj <- readRDS(file_path)
    
    # Validate object
    validation_errors <- validate_seurat_object(sobj, file_path)
    
    if (length(validation_errors) > 0) {
      cat("Skipping", lineage_name, "due to validation errors\n")
      next
    }
    
    cat("Loaded:", ncol(sobj), "cells,", nrow(sobj), "genes\n")
    
    # Perform comparative analysis for all condition pairs
    if (length(CONDITIONS) >= 2) {
      cond1 <- CONDITIONS[1]
      cond2 <- CONDITIONS[2]
      
      deg_results <- perform_deg_analysis(sobj, cond1, cond2, lineage_name)
      
      if (!is.null(deg_results)) {
        all_lineage_results[[lineage_name]] <- deg_results
        cat("Successfully analyzed", lineage_name, "\n")
      } else {
        cat("No results for", lineage_name, "\n")
      }
    }
    
  }, error = function(e) {
    cat("Error processing", lineage_name, ":", e$message, "\n")
  })
}

# ============================================================================
# Combine and Analyze Results
# ============================================================================

cat("\nCombining results from all lineages...\n")

if (length(all_lineage_results) > 0) {
  # Combine all results
  all_results <- bind_rows(all_lineage_results)
  cat("Combined results:", nrow(all_results), "DEGs from", length(all_lineage_results), "lineages\n")
  
  # Create summary
  summary_data <- create_comparison_summary(all_results)
  
  # Print summary
  if (!is.null(summary_data$overall)) {
    cat("\nOverall Summary:\n")
    print(summary_data$overall)
  }
  
  if (!is.null(summary_data$lineage)) {
    cat("\nTop 5 Lineages by DEG Count:\n")
    print(head(summary_data$lineage, 5))
  }
  
} else {
  cat("No results found from any lineage\n")
  all_results <- NULL
  summary_data <- NULL
}

# ============================================================================
# Create Visualizations and Save Results
# ============================================================================

if (!is.null(all_results) && nrow(all_results) > 0) {
  
  cat("\nCreating visualizations...\n")
  
  # Define output prefix
  output_prefix <- paste0("multi_comparison_", paste(CONDITIONS, collapse = "_vs_"), "_", JOB_ID)
  
  # Create cross-lineage plots
  plots_created <- create_cross_lineage_plots(all_results, output_prefix)
  
  # Save comprehensive results
  files_saved <- save_comprehensive_results(all_results, summary_data, output_prefix)
  
  cat("Files saved:", length(files_saved), "\n")
  cat("Plots created:", length(plots_created), "\n")
  
} else {
  cat("No results to save or visualize\n")
  files_saved <- character(0)
  plots_created <- character(0)
}

# ============================================================================
# Generate Report
# ============================================================================

cat("\nGenerating comprehensive report...\n")

report_text <- paste0(
  "Multi-Condition Comparative Analysis Report\n",
  "==========================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Job ID: ", JOB_ID, "\n",
  "Conditions Compared: ", paste(CONDITIONS, collapse = " vs "), "\n",
  "Input Directory: ", INPUT_DIR, "\n",
  "Files Processed: ", length(input_files), "\n",
  "Successful Lineages: ", length(all_lineage_results), "\n\n"
)

if (!is.null(summary_data$overall)) {
  report_text <- paste0(report_text,
    "Overall Results:\n",
    "  - Total lineages analyzed: ", summary_data$overall$total_lineages, "\n",
    "  - Total clusters: ", summary_data$overall$total_clusters, "\n",
    "  - Total DEGs found: ", summary_data$overall$total_genes, "\n",
    "  - Significant DEGs: ", summary_data$overall$significant_genes, "\n",
    "  - Mean |Log2FC|: ", round(summary_data$overall$mean_log2FC, 3), "\n\n"
  )
}

if (!is.null(summary_data$lineage) && nrow(summary_data$lineage) > 0) {
  report_text <- paste0(report_text, "Top 5 Lineages by DEG Count:\n")
  top_lineages <- head(summary_data$lineage, 5)
  for (i in 1:nrow(top_lineages)) {
    lineage_info <- top_lineages[i, ]
    report_text <- paste0(report_text,
      "  ", i, ". ", lineage_info$lineage, ": ", lineage_info$total_genes, 
      " DEGs (", lineage_info$significant_genes, " significant)\n"
    )
  }
  report_text <- paste0(report_text, "\n")
}

report_text <- paste0(report_text,
  "Files Generated:\n"
)

if (length(files_saved) > 0) {
  for (file in files_saved) {
    report_text <- paste0(report_text, "  - ", basename(file), "\n")
  }
}

report_text <- paste0(report_text, "\nPlots Generated:\n")
if (length(plots_created) > 0) {
  for (plot in plots_created) {
    report_text <- paste0(report_text, "  - ", plot, ".pdf\n")
  }
}

report_text <- paste0(report_text,
  "\nRecommendations:\n",
  "1. Focus on lineages with highest DEG counts for validation\n",
  "2. Investigate common genes across multiple lineages\n", 
  "3. Perform functional enrichment analysis on DEG sets\n",
  "4. Consider batch effects if results vary unexpectedly\n",
  "5. Validate key findings with independent methods\n\n",
  "Usage Examples:\n",
  "Rscript multi_condition_comparison.R <input_dir> <job_id> <cond1,cond2>\n",
  "Rscript multi_condition_comparison.R data/ analysis_v1 CS9,CS10\n"
)

# Save report
report_file <- file.path(OUTPUT_DIR, paste0("multi_comparison_report_", JOB_ID, ".txt"))
writeLines(report_text, report_file)

# Save session info
session_file <- file.path(OUTPUT_DIR, paste0("multi_comparison_sessionInfo_", JOB_ID, ".txt"))
writeLines(capture.output(sessionInfo()), session_file)

# ============================================================================
# Final Summary
# ============================================================================

cat("\n", "="*60, "\n")
cat("Multi-condition comparative analysis completed!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\nFinal Summary:\n")
cat("- Input files processed:", length(input_files), "\n")
cat("- Successful lineages:", length(all_lineage_results), "\n")
cat("- Conditions compared:", paste(CONDITIONS, collapse = " vs "), "\n")

if (!is.null(all_results)) {
  cat("- Total DEGs found:", nrow(all_results), "\n")
  cat("- Lineages with DEGs:", n_distinct(all_results$lineage), "\n")
}

cat("- Files saved:", length(files_saved), "\n")
cat("- Plots created:", length(plots_created), "\n")

cat("\nNext steps:\n")
cat("1. Review lineage-specific results for biological insights\n")
cat("2. Identify conserved and lineage-specific DEG patterns\n")
cat("3. Perform cross-lineage functional analysis\n")
cat("4. Validate findings with experimental approaches\n")

cat("\nAnalysis completed successfully!\n")