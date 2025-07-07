#!/usr/bin/env Rscript

# ============================================================================
# Comprehensive Differential Expression Analysis Between Conditions
# ============================================================================
# 
# This script performs systematic differential expression analysis comparing
# different conditions, developmental stages, or treatments across multiple
# cell types and clusters in spatial transcriptomics data.
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
  library(EnhancedVolcano)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(future.apply)
})

# Parameters
JOB_ID <- "deg_analysis"
CONDITION_1 <- "CS9"              # First condition to compare
CONDITION_2 <- "CS10"             # Second condition to compare
CONDITION_COL <- "CS"             # Column name for conditions
CELL_TYPE_COL <- "celltype"       # Column name for cell types
MIN_PCT <- 0.25                   # Minimum percentage of cells expressing gene
LOGFC_THRESHOLD <- 0.25           # Log fold change threshold
ONLY_POS <- FALSE                 # Include both up and down regulated genes
P_VAL_CUTOFF <- 0.05              # P-value cutoff for significance
TOP_N_GENES <- 100                # Number of top genes to highlight
N_WORKERS <- 4                    # Number of parallel workers

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
  CONDITION_1 <- args[2]
}
if (length(args) > 2) {
  CONDITION_2 <- args[3]
}
if (length(args) > 3) {
  JOB_ID <- args[4]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to perform DEG analysis for a specific cluster
get_markers_for_cluster <- function(cluster_id, sobj, cluster_sizes, 
                                   condition_col = CONDITION_COL,
                                   cond1 = CONDITION_1, cond2 = CONDITION_2) {
  
  # Check if cluster has enough cells
  if (cluster_sizes[[cluster_id]] < 3) {
    message("Skipping cluster ", cluster_id, ": fewer than 3 cells")
    return(NULL)
  }
  
  # Subset to specific cluster
  sobj_cluster <- subset(sobj, idents = cluster_id)
  
  # Check if both conditions are present
  if (!all(c(cond1, cond2) %in% unique(sobj_cluster@meta.data[[condition_col]]))) {
    message("Skipping cluster ", cluster_id, ": missing ", cond1, " or ", cond2)
    return(NULL)
  }
  
  # Check cell counts per condition
  condition_counts <- table(sobj_cluster@meta.data[[condition_col]])
  if (any(condition_counts < 3)) {
    message("Skipping cluster ", cluster_id, ": too few cells in one condition")
    return(NULL)
  }
  
  # Perform differential expression analysis
  tryCatch({
    markers <- FindMarkers(
      sobj_cluster, 
      ident.1 = cond1, 
      ident.2 = cond2, 
      group.by = condition_col, 
      only.pos = ONLY_POS, 
      min.pct = MIN_PCT, 
      logfc.threshold = LOGFC_THRESHOLD,
      test.use = "wilcox"
    )
    
    # Process results
    if (nrow(markers) > 0) {
      markers <- markers %>%
        as.data.frame() %>%
        mutate(
          gene = rownames(.),
          cluster = cluster_id,
          comparison = paste(cond1, "vs", cond2),
          significant = ifelse(p_val_adj < P_VAL_CUTOFF, "Significant", "Not significant"),
          regulation = ifelse(avg_log2FC > 0, paste("Up in", cond1), paste("Up in", cond2)),
          abs_log2FC = abs(avg_log2FC),
          neg_log10_pval = -log10(p_val_adj + 1e-300)
        ) %>%
        select(gene, cluster, comparison, everything()) %>%
        arrange(desc(abs_log2FC))
      
      message("Found ", nrow(markers), " DEGs for cluster ", cluster_id)
      return(markers)
    } else {
      message("No DEGs found for cluster ", cluster_id)
      return(NULL)
    }
    
  }, error = function(e) {
    message("Error analyzing cluster ", cluster_id, ": ", e$message)
    return(NULL)
  })
}

# Function to create Excel file with multiple sheets
create_deg_excel <- function(markers_list, filename, zoom = 200) {
  wb <- createWorkbook()
  
  # Create summary sheet
  if (length(markers_list) > 0) {
    # Combine all markers for summary
    all_markers <- bind_rows(markers_list)
    
    # Summary statistics
    summary_stats <- all_markers %>%
      group_by(cluster) %>%
      summarise(
        total_genes = n(),
        up_regulated = sum(regulation == paste("Up in", CONDITION_1)),
        down_regulated = sum(regulation == paste("Up in", CONDITION_2)),
        significant = sum(significant == "Significant"),
        mean_log2FC = mean(abs_log2FC),
        .groups = "drop"
      ) %>%
      arrange(desc(total_genes))
    
    # Add summary sheet
    addWorksheet(wb, "Summary", zoom = zoom)
    writeData(wb, "Summary", summary_stats)
    
    # Add all markers sheet
    addWorksheet(wb, "All_Markers", zoom = zoom)
    writeData(wb, "All_Markers", all_markers)
  }
  
  # Add individual cluster sheets
  sheet_names <- substr(paste0("Cluster_", names(markers_list)), 1, 31)
  
  for (i in seq_along(markers_list)) {
    if (!is.null(markers_list[[i]]) && nrow(markers_list[[i]]) > 0) {
      addWorksheet(wb, sheet_names[i], zoom = zoom)
      writeData(wb, sheet_names[i], markers_list[[i]])
    }
  }
  
  saveWorkbook(wb, filename, overwrite = TRUE)
  message("Excel file saved: ", filename)
}

# Function to create volcano plot
create_volcano_plot <- function(markers, cluster_name, output_file) {
  
  if (is.null(markers) || nrow(markers) == 0) {
    message("No data for volcano plot: ", cluster_name)
    return(NULL)
  }
  
  # Prepare data for volcano plot
  plot_data <- markers %>%
    mutate(
      log2FC = avg_log2FC,
      neg_log10_pval = -log10(p_val_adj + 1e-300),
      label = ifelse(abs(log2FC) > 1 & p_val_adj < 0.01, gene, "")
    )
  
  # Create volcano plot
  p <- ggplot(plot_data, aes(x = log2FC, y = neg_log10_pval)) +
    geom_point(aes(color = regulation), alpha = 0.6, size = 1) +
    geom_hline(yintercept = -log10(P_VAL_CUTOFF), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Up in CS9" = "red", "Up in CS10" = "blue")) +
    labs(
      title = paste("Volcano Plot -", cluster_name),
      subtitle = paste(CONDITION_1, "vs", CONDITION_2),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Regulation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Add gene labels for top genes
  if (any(plot_data$label != "")) {
    p <- p + geom_text_repel(aes(label = label), size = 3, max.overlaps = 20)
  }
  
  # Save plot
  ggsave(output_file, p, width = 10, height = 8, dpi = 300)
  message("Volcano plot saved: ", output_file)
  
  return(p)
}

# Function to create heatmap of top DEGs
create_deg_heatmap <- function(sobj, all_markers, output_file, top_n = 50) {
  
  if (is.null(all_markers) || nrow(all_markers) == 0) {
    message("No data for heatmap")
    return(NULL)
  }
  
  # Select top genes
  top_genes <- all_markers %>%
    arrange(desc(abs_log2FC)) %>%
    head(top_n) %>%
    pull(gene)
  
  # Aggregate expression by condition and cell type
  sobj$condition_celltype <- paste(sobj@meta.data[[CONDITION_COL]], 
                                   sobj@meta.data[[CELL_TYPE_COL]], 
                                   sep = "_")
  
  # Calculate average expression
  avg_exp <- AverageExpression(sobj, 
                              features = top_genes, 
                              group.by = "condition_celltype",
                              assays = "RNA")
  
  # Extract expression matrix
  exp_matrix <- as.matrix(avg_exp$RNA)
  
  # Create heatmap
  pdf(output_file, width = 12, height = 10)
  
  pheatmap(exp_matrix,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           main = paste("Top", top_n, "DEGs Heatmap"),
           fontsize_row = 8,
           fontsize_col = 10,
           angle_col = 45)
  
  dev.off()
  
  message("DEG heatmap saved: ", output_file)
  
  return(exp_matrix)
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

cat("Starting differential expression analysis...\n")
cat("Comparison:", CONDITION_1, "vs", CONDITION_2, "\n")
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
  } else if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
    CELL_TYPE_COL <- "seurat_clusters"
  } else {
    stop("Cell type column not found")
  }
}

cat("Using condition column:", CONDITION_COL, "\n")
cat("Using cell type column:", CELL_TYPE_COL, "\n")

# Check conditions
available_conditions <- unique(sobj@meta.data[[CONDITION_COL]])
cat("Available conditions:", paste(available_conditions, collapse = ", "), "\n")

if (!CONDITION_1 %in% available_conditions) {
  stop("Condition 1 not found: ", CONDITION_1)
}
if (!CONDITION_2 %in% available_conditions) {
  stop("Condition 2 not found: ", CONDITION_2)
}

# Set up parallel processing
plan(multisession, workers = N_WORKERS)
options(future.globals.maxSize = 3 * 1024^3)

# Set cell type as identity
Idents(sobj) <- CELL_TYPE_COL

# Get clusters and their sizes
clusters <- Idents(sobj) %>% unique() %>% sort()
cluster_sizes <- table(Idents(sobj))

cat("Found", length(clusters), "cell types/clusters\n")
cat("Cell type distribution:\n")
print(cluster_sizes)

# ============================================================================
# Perform DEG Analysis
# ============================================================================

cat("Performing DEG analysis for each cluster...\n")

# Perform analysis for each cluster
markers_list <- clusters %>%
  set_names() %>%
  map(~ get_markers_for_cluster(.x, sobj, cluster_sizes)) %>%
  discard(is.null)

cat("Successfully analyzed", length(markers_list), "clusters\n")

# Check if any results were found
if (length(markers_list) == 0) {
  stop("No DEGs found for any cluster. Check your conditions and parameters.")
}

# Combine all markers
all_markers <- bind_rows(markers_list)
cat("Total DEGs found:", nrow(all_markers), "\n")

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\n")

# Define output file names
output_prefix <- paste0("DEG_", CONDITION_1, "_vs_", CONDITION_2, "_", JOB_ID)

# Save individual results
rds_file <- file.path(DATA_DIR, paste0(output_prefix, ".rds"))
csv_file <- file.path(DATA_DIR, paste0(output_prefix, ".csv"))
xlsx_file <- file.path(OUTPUT_DIR, paste0(output_prefix, ".xlsx"))

# Save as RDS and CSV
saveRDS(all_markers, rds_file)
write.csv(all_markers, csv_file, row.names = FALSE)

# Create comprehensive Excel file
create_deg_excel(markers_list, xlsx_file)

# ============================================================================
# Create Visualizations
# ============================================================================

cat("Creating visualizations...\n")

# 1. Create summary statistics plot
summary_stats <- all_markers %>%
  group_by(cluster) %>%
  summarise(
    total_genes = n(),
    up_regulated = sum(regulation == paste("Up in", CONDITION_1)),
    down_regulated = sum(regulation == paste("Up in", CONDITION_2)),
    significant = sum(significant == "Significant"),
    mean_log2FC = mean(abs_log2FC),
    .groups = "drop"
  ) %>%
  arrange(desc(total_genes))

# Summary barplot
p_summary <- summary_stats %>%
  pivot_longer(cols = c(up_regulated, down_regulated), 
               names_to = "direction", values_to = "count") %>%
  ggplot(aes(x = reorder(cluster, total_genes), y = count, fill = direction)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c("up_regulated" = "red", "down_regulated" = "blue"),
                    labels = c(paste("Up in", CONDITION_1), paste("Up in", CONDITION_2))) +
  labs(title = "DEG Summary by Cluster",
       x = "Cluster", y = "Number of DEGs", fill = "Direction") +
  theme_minimal()

ggsave(file.path(PLOTS_DIR, paste0(output_prefix, "_summary.pdf")), 
       p_summary, width = 10, height = 8, dpi = 300)

# 2. Create volcano plots for top clusters
top_clusters <- head(summary_stats$cluster, 5)

for (cluster in top_clusters) {
  volcano_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_volcano_", cluster, ".pdf"))
  create_volcano_plot(markers_list[[cluster]], cluster, volcano_file)
}

# 3. Create overall heatmap
heatmap_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_heatmap.pdf"))
create_deg_heatmap(sobj, all_markers, heatmap_file, top_n = TOP_N_GENES)

# 4. Create pathway enrichment summary (if available)
if (nrow(all_markers) > 0) {
  # Top genes by fold change
  top_genes_plot <- all_markers %>%
    arrange(desc(abs_log2FC)) %>%
    head(20) %>%
    ggplot(aes(x = reorder(gene, abs_log2FC), y = abs_log2FC, fill = regulation)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("red", "blue")) +
    labs(title = "Top 20 DEGs by Fold Change",
         x = "Gene", y = "Absolute Log2 Fold Change") +
    theme_minimal()
  
  ggsave(file.path(PLOTS_DIR, paste0(output_prefix, "_top_genes.pdf")), 
         top_genes_plot, width = 10, height = 8, dpi = 300)
}

# ============================================================================
# Generate Report
# ============================================================================

cat("Generating analysis report...\n")

# Calculate overall statistics
n_total_genes <- nrow(all_markers)
n_significant <- sum(all_markers$significant == "Significant")
n_up_cond1 <- sum(all_markers$regulation == paste("Up in", CONDITION_1))
n_up_cond2 <- sum(all_markers$regulation == paste("Up in", CONDITION_2))
n_clusters_analyzed <- length(markers_list)

# Generate report
report_text <- paste0(
  "Differential Expression Analysis Report\n",
  "=====================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Job ID: ", JOB_ID, "\n",
  "Comparison: ", CONDITION_1, " vs ", CONDITION_2, "\n\n",
  "Parameters:\n",
  "  - Minimum percentage: ", MIN_PCT, "\n",
  "  - LogFC threshold: ", LOGFC_THRESHOLD, "\n",
  "  - P-value cutoff: ", P_VAL_CUTOFF, "\n",
  "  - Only positive: ", ONLY_POS, "\n\n",
  "Results Summary:\n",
  "  - Total cells analyzed: ", ncol(sobj), "\n",
  "  - Clusters analyzed: ", n_clusters_analyzed, "\n",
  "  - Total DEGs found: ", n_total_genes, "\n",
  "  - Significant DEGs: ", n_significant, " (", round(n_significant/n_total_genes*100, 1), "%)\n",
  "  - Up in ", CONDITION_1, ": ", n_up_cond1, "\n",
  "  - Up in ", CONDITION_2, ": ", n_up_cond2, "\n\n",
  "Top 5 Clusters by DEG Count:\n"
)

for (i in 1:min(5, nrow(summary_stats))) {
  cluster_info <- summary_stats[i, ]
  report_text <- paste0(report_text,
    "  ", i, ". ", cluster_info$cluster, ": ", cluster_info$total_genes, " DEGs\n"
  )
}

report_text <- paste0(report_text,
  "\nTop 10 DEGs by Fold Change:\n"
)

top_10_degs <- all_markers %>% arrange(desc(abs_log2FC)) %>% head(10)
for (i in 1:min(10, nrow(top_10_degs))) {
  deg_info <- top_10_degs[i, ]
  report_text <- paste0(report_text,
    "  ", i, ". ", deg_info$gene, " (", deg_info$cluster, "): ", 
    round(deg_info$avg_log2FC, 3), " log2FC\n"
  )
}

report_text <- paste0(report_text,
  "\nOutput Files:\n",
  "  - Combined results: ", basename(csv_file), "\n",
  "  - Excel workbook: ", basename(xlsx_file), "\n",
  "  - R object: ", basename(rds_file), "\n\n",
  "Visualizations:\n",
  "  - Summary plot: ", output_prefix, "_summary.pdf\n",
  "  - Top genes plot: ", output_prefix, "_top_genes.pdf\n",
  "  - DEG heatmap: ", output_prefix, "_heatmap.pdf\n",
  "  - Volcano plots: ", output_prefix, "_volcano_[cluster].pdf\n\n",
  "Recommendations:\n",
  "1. Validate top DEGs with qPCR or other methods\n",
  "2. Perform functional enrichment analysis (GO, KEGG)\n",
  "3. Check for batch effects if results seem unexpected\n",
  "4. Consider pathway analysis for biological interpretation\n",
  "5. Examine spatial distribution of DEGs if applicable\n\n",
  "Usage Examples:\n",
  "Rscript differential_expression_analysis.R <input.rds> <cond1> <cond2> <job_id>\n",
  "Rscript differential_expression_analysis.R data.rds CS9 CS10 analysis_v1\n"
)

# Save report
report_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_report.txt"))
writeLines(report_text, report_file)

# Save session info
session_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), session_file)

# ============================================================================
# Final Summary
# ============================================================================

cat("\n", "="*60, "\n")
cat("Differential expression analysis completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")
cat("\nFinal Summary:\n")
cat("- Clusters analyzed:", n_clusters_analyzed, "\n")
cat("- Total DEGs found:", n_total_genes, "\n")
cat("- Significant DEGs:", n_significant, "\n")
cat("- Up-regulated in", CONDITION_1, ":", n_up_cond1, "\n")
cat("- Up-regulated in", CONDITION_2, ":", n_up_cond2, "\n")

if (n_total_genes > 0) {
  top_cluster <- summary_stats$cluster[1]
  cat("- Top cluster by DEG count:", top_cluster, 
      "(", summary_stats$total_genes[1], "DEGs )\n")
}

cat("\nNext steps:\n")
cat("1. Review DEG results and volcano plots\n")
cat("2. Perform functional enrichment analysis\n")
cat("3. Validate key findings experimentally\n")
cat("4. Compare with literature and databases\n")

cat("\nAnalysis completed successfully!\n")