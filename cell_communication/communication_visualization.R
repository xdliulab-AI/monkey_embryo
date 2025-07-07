#!/usr/bin/env Rscript

# ============================================================================
# Cell-Cell Communication Visualization and Comparison
# ============================================================================
# 
# This script creates comprehensive visualizations for cell-cell communication
# analysis results from CellChat and other methods. It includes dot plots,
# bubble plots, network visualizations, and comparative analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(scales)
  library(gridExtra)
  library(ComplexHeatmap)
  library(circlize)
  library(pheatmap)
  library(corrplot)
})

# Parameters
INPUT_CSV <- "../data/ligand_receptor_interactions.csv"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"
JOB_ID <- "communication_viz"

# Visualization parameters
POINT_SIZE_RANGE <- c(1, 8)
ALPHA_VALUE <- 0.7
TEXT_ANGLE <- 45
BASE_TEXT_SIZE <- 12
DPI_VALUE <- 300
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 10

# Color palettes
PALETTES <- list(
  viridis = viridis_d(10),
  plasma = plasma(10),
  inferno = inferno(10),
  magma = magma(10),
  blues = brewer.pal(9, "Blues"),
  reds = brewer.pal(9, "Reds"),
  greens = brewer.pal(9, "Greens"),
  spectral = brewer.pal(11, "Spectral"),
  rdylbu = brewer.pal(11, "RdYlBu"),
  set3 = brewer.pal(12, "Set3")
)

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  INPUT_CSV <- args[1]
}
if (length(args) > 1) {
  JOB_ID <- args[2]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to clean p-values
clean_pvalues <- function(pval_vector, min_val = 1e-10) {
  pval_vector[pval_vector == 0] <- min_val
  pval_vector[is.na(pval_vector)] <- 1
  return(pval_vector)
}

# Function to create source-target labels
create_interaction_labels <- function(df) {
  df$source_target <- paste(df$source, df$target, sep = " → ")
  df$ligand_receptor <- paste(df$ligand, df$receptor, sep = "-")
  
  # Use interaction_name_2 if available, otherwise create from ligand-receptor
  if ("interaction_name_2" %in% colnames(df)) {
    df$lr_pair <- df$interaction_name_2
  } else if ("interaction_name" %in% colnames(df)) {
    df$lr_pair <- df$interaction_name
  } else {
    df$lr_pair <- df$ligand_receptor
  }
  
  return(df)
}

# Function to calculate significance levels
add_significance_levels <- function(df, pval_col = "pval") {
  df$significance <- case_when(
    df[[pval_col]] < 0.001 ~ "***",
    df[[pval_col]] < 0.01 ~ "**",
    df[[pval_col]] < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  df$neg_log_pval <- -log10(df[[pval_col]])
  return(df)
}

# Function to filter top interactions
filter_top_interactions <- function(df, top_n = 50, by_column = "prob") {
  if (nrow(df) > top_n) {
    df <- df %>%
      arrange(desc(.data[[by_column]])) %>%
      slice_head(n = top_n)
    cat("Filtered to top", top_n, "interactions by", by_column, "\n")
  }
  return(df)
}

# Function to create basic dot plot
create_dotplot <- function(df, palette_name = "viridis", 
                          size_col = "neg_log_pval", color_col = "prob",
                          max_interactions = 100) {
  
  # Filter if too many interactions
  if (nrow(df) > max_interactions) {
    df <- filter_top_interactions(df, max_interactions, color_col)
  }
  
  # Convert to factors for proper ordering
  df$source_target <- factor(df$source_target)
  df$lr_pair <- factor(df$lr_pair)
  
  # Create base plot
  p <- ggplot(df, aes(x = source_target, y = lr_pair)) +
    geom_point(aes_string(color = color_col, size = size_col), alpha = ALPHA_VALUE) +
    scale_size_continuous(
      range = POINT_SIZE_RANGE,
      name = expression(-log[10](p)),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    theme_bw(base_size = BASE_TEXT_SIZE) +
    theme(
      axis.text.x = element_text(angle = TEXT_ANGLE, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      x = "Cell Type Pair (Source → Target)",
      y = "Ligand-Receptor Pair",
      title = paste("Cell-Cell Communication:", palette_name, "palette")
    )
  
  # Add color scale based on palette
  if (palette_name %in% names(PALETTES)) {
    p <- p + scale_color_gradientn(
      colors = PALETTES[[palette_name]],
      name = "Comm. Prob.",
      guide = guide_colorbar(override.aes = list(alpha = 1))
    )
  } else {
    p <- p + scale_color_viridis_c(
      name = "Comm. Prob.",
      option = "viridis"
    )
  }
  
  return(p)
}

# Function to create bubble plot
create_bubble_plot <- function(df, max_interactions = 80) {
  
  # Filter if too many interactions
  if (nrow(df) > max_interactions) {
    df <- filter_top_interactions(df, max_interactions, "prob")
  }
  
  # Convert to factors
  df$source_target <- factor(df$source_target)
  df$lr_pair <- factor(df$lr_pair)
  
  p <- ggplot(df, aes(x = source_target, y = lr_pair)) +
    geom_point(aes(size = neg_log_pval, color = prob), alpha = ALPHA_VALUE) +
    scale_size_continuous(
      name = expression(-log[10](p)),
      range = c(1, 10),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_color_viridis_c(
      name = "Probability",
      option = "viridis"
    ) +
    theme_minimal(base_size = BASE_TEXT_SIZE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "gray", size = 0.2),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Ligand-Receptor Interaction Bubble Plot",
      x = "Source - Target Cell Type",
      y = "Interaction (Ligand - Receptor)"
    )
  
  return(p)
}

# Function to create communication heatmap by cell type
create_celltype_heatmap <- function(df) {
  
  # Aggregate by cell type pairs
  comm_summary <- df %>%
    group_by(source, target) %>%
    summarise(
      mean_prob = mean(prob, na.rm = TRUE),
      total_interactions = n(),
      mean_significance = mean(neg_log_pval, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create matrix for heatmap
  cell_types <- unique(c(comm_summary$source, comm_summary$target))
  comm_matrix <- matrix(0, nrow = length(cell_types), ncol = length(cell_types))
  rownames(comm_matrix) <- cell_types
  colnames(comm_matrix) <- cell_types
  
  for (i in 1:nrow(comm_summary)) {
    row_idx <- which(cell_types == comm_summary$source[i])
    col_idx <- which(cell_types == comm_summary$target[i])
    comm_matrix[row_idx, col_idx] <- comm_summary$mean_prob[i]
  }
  
  # Create heatmap
  p <- pheatmap(
    comm_matrix,
    color = viridis(100),
    display_numbers = TRUE,
    number_format = "%.3f",
    fontsize_number = 8,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Mean Communication Probability by Cell Type",
    silent = TRUE
  )
  
  return(list(plot = p, matrix = comm_matrix, summary = comm_summary))
}

# Function to create pathway analysis
create_pathway_analysis <- function(df) {
  
  if (!"pathway_name" %in% colnames(df)) {
    cat("Warning: No pathway information found\n")
    return(NULL)
  }
  
  # Summarize by pathway
  pathway_summary <- df %>%
    group_by(pathway_name) %>%
    summarise(
      n_interactions = n(),
      mean_prob = mean(prob, na.rm = TRUE),
      max_prob = max(prob, na.rm = TRUE),
      mean_significance = mean(neg_log_pval, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_prob))
  
  # Create pathway plot
  top_pathways <- head(pathway_summary, 20)
  
  p <- ggplot(top_pathways, aes(x = reorder(pathway_name, mean_prob), y = mean_prob)) +
    geom_col(aes(fill = mean_significance), alpha = 0.8) +
    scale_fill_viridis_c(name = "Mean\n-log10(p)") +
    coord_flip() +
    theme_minimal(base_size = BASE_TEXT_SIZE) +
    labs(
      title = "Top Communication Pathways",
      x = "Pathway",
      y = "Mean Communication Probability"
    ) +
    theme(
      axis.text.y = element_text(size = 10)
    )
  
  return(list(plot = p, summary = pathway_summary))
}

# Function to create network summary plot
create_network_summary <- function(df) {
  
  # Calculate outgoing and incoming communication per cell type
  outgoing <- df %>%
    group_by(source) %>%
    summarise(
      outgoing_strength = sum(prob, na.rm = TRUE),
      outgoing_count = n(),
      .groups = "drop"
    ) %>%
    rename(cell_type = source)
  
  incoming <- df %>%
    group_by(target) %>%
    summarise(
      incoming_strength = sum(prob, na.rm = TRUE),
      incoming_count = n(),
      .groups = "drop"
    ) %>%
    rename(cell_type = target)
  
  # Merge outgoing and incoming
  network_stats <- full_join(outgoing, incoming, by = "cell_type") %>%
    mutate(
      outgoing_strength = ifelse(is.na(outgoing_strength), 0, outgoing_strength),
      incoming_strength = ifelse(is.na(incoming_strength), 0, incoming_strength),
      outgoing_count = ifelse(is.na(outgoing_count), 0, outgoing_count),
      incoming_count = ifelse(is.na(incoming_count), 0, incoming_count),
      total_strength = outgoing_strength + incoming_strength,
      net_strength = outgoing_strength - incoming_strength
    )
  
  # Create visualization
  p1 <- ggplot(network_stats, aes(x = reorder(cell_type, total_strength))) +
    geom_col(aes(y = outgoing_strength), fill = "red", alpha = 0.7) +
    geom_col(aes(y = -incoming_strength), fill = "blue", alpha = 0.7) +
    coord_flip() +
    theme_minimal(base_size = BASE_TEXT_SIZE) +
    labs(
      title = "Communication Strength by Cell Type",
      subtitle = "Red: Outgoing, Blue: Incoming",
      x = "Cell Type",
      y = "Communication Strength"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(list(plot = p1, stats = network_stats))
}

# ============================================================================
# Main Analysis
# ============================================================================

cat("Loading and processing communication data...\n")

# Load data
if (!file.exists(INPUT_CSV)) {
  stop("Input CSV file not found: ", INPUT_CSV)
}

df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
cat("Loaded data with", nrow(df), "interactions\n")

# Display column names
cat("Available columns:", paste(colnames(df), collapse = ", "), "\n")

# Process data
df <- create_interaction_labels(df)
df$pval <- clean_pvalues(df$pval)
df <- add_significance_levels(df)

cat("Data processing completed\n")
cat("Unique source cell types:", length(unique(df$source)), "\n")
cat("Unique target cell types:", length(unique(df$target)), "\n")
cat("Unique L-R pairs:", length(unique(df$lr_pair)), "\n")

# ============================================================================
# Create Visualizations
# ============================================================================

cat("Creating visualizations...\n")

# 1. Create dot plots with different palettes
cat("Creating dot plots with different color palettes...\n")

for (palette_name in names(PALETTES)) {
  tryCatch({
    p <- create_dotplot(df, palette_name = palette_name, max_interactions = 80)
    
    filename <- file.path(PLOTS_DIR, 
                         paste0("dotplot_", palette_name, "_", JOB_ID, ".pdf"))
    
    ggsave(filename, p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI_VALUE)
    cat("Saved:", basename(filename), "\n")
    
  }, error = function(e) {
    cat("Error creating", palette_name, "dot plot:", e$message, "\n")
  })
}

# 2. Create bubble plot
cat("Creating bubble plot...\n")
tryCatch({
  bubble_plot <- create_bubble_plot(df, max_interactions = 100)
  
  filename <- file.path(PLOTS_DIR, paste0("bubble_plot_", JOB_ID, ".pdf"))
  ggsave(filename, bubble_plot, width = PLOT_WIDTH, height = 14, dpi = DPI_VALUE)
  cat("Saved:", basename(filename), "\n")
  
}, error = function(e) {
  cat("Error creating bubble plot:", e$message, "\n")
})

# 3. Create cell type communication heatmap
cat("Creating cell type communication heatmap...\n")
tryCatch({
  heatmap_result <- create_celltype_heatmap(df)
  
  filename <- file.path(PLOTS_DIR, paste0("celltype_heatmap_", JOB_ID, ".pdf"))
  pdf(filename, width = 10, height = 8)
  print(heatmap_result$plot)
  dev.off()
  cat("Saved:", basename(filename), "\n")
  
  # Save communication matrix
  write.csv(heatmap_result$matrix, 
            file.path(OUTPUT_DIR, paste0("communication_matrix_", JOB_ID, ".csv")))
  write.csv(heatmap_result$summary, 
            file.path(OUTPUT_DIR, paste0("celltype_communication_summary_", JOB_ID, ".csv")),
            row.names = FALSE)
  
}, error = function(e) {
  cat("Error creating heatmap:", e$message, "\n")
})

# 4. Create pathway analysis
cat("Creating pathway analysis...\n")
tryCatch({
  pathway_result <- create_pathway_analysis(df)
  
  if (!is.null(pathway_result)) {
    filename <- file.path(PLOTS_DIR, paste0("pathway_analysis_", JOB_ID, ".pdf"))
    ggsave(filename, pathway_result$plot, width = 10, height = 8, dpi = DPI_VALUE)
    cat("Saved:", basename(filename), "\n")
    
    # Save pathway summary
    write.csv(pathway_result$summary, 
              file.path(OUTPUT_DIR, paste0("pathway_summary_", JOB_ID, ".csv")),
              row.names = FALSE)
  }
  
}, error = function(e) {
  cat("Error creating pathway analysis:", e$message, "\n")
})

# 5. Create network summary
cat("Creating network summary...\n")
tryCatch({
  network_result <- create_network_summary(df)
  
  filename <- file.path(PLOTS_DIR, paste0("network_summary_", JOB_ID, ".pdf"))
  ggsave(filename, network_result$plot, width = 10, height = 8, dpi = DPI_VALUE)
  cat("Saved:", basename(filename), "\n")
  
  # Save network statistics
  write.csv(network_result$stats, 
            file.path(OUTPUT_DIR, paste0("network_statistics_", JOB_ID, ".csv")),
            row.names = FALSE)
  
}, error = function(e) {
  cat("Error creating network summary:", e$message, "\n")
})

# ============================================================================
# Generate Summary Statistics
# ============================================================================

cat("Generating summary statistics...\n")

# Calculate overall statistics
total_interactions <- nrow(df)
unique_lr_pairs <- length(unique(df$lr_pair))
unique_cell_pairs <- length(unique(df$source_target))
mean_prob <- mean(df$prob, na.rm = TRUE)
median_prob <- median(df$prob, na.rm = TRUE)

# Identify top interactions
top_interactions <- df %>%
  arrange(desc(prob)) %>%
  head(10) %>%
  select(source, target, lr_pair, prob, pval, significance)

# Identify top cell type pairs
top_cell_pairs <- df %>%
  group_by(source, target) %>%
  summarise(
    mean_prob = mean(prob, na.rm = TRUE),
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_prob)) %>%
  head(10)

# ============================================================================
# Generate Report
# ============================================================================

cat("Generating comprehensive report...\n")

report_text <- paste0(
  "Cell-Cell Communication Visualization Analysis Report\n",
  "===================================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Job ID: ", JOB_ID, "\n",
  "Input File: ", basename(INPUT_CSV), "\n\n",
  "Data Summary:\n",
  "  - Total interactions: ", total_interactions, "\n",
  "  - Unique L-R pairs: ", unique_lr_pairs, "\n",
  "  - Unique cell type pairs: ", unique_cell_pairs, "\n",
  "  - Mean communication probability: ", round(mean_prob, 4), "\n",
  "  - Median communication probability: ", round(median_prob, 4), "\n\n",
  "Top 5 Interactions by Probability:\n"
)

for (i in 1:min(5, nrow(top_interactions))) {
  interaction <- top_interactions[i, ]
  report_text <- paste0(report_text,
    "  ", i, ". ", interaction$source, " → ", interaction$target, 
    ": ", interaction$lr_pair, " (prob: ", round(interaction$prob, 4), 
    ", p: ", format(interaction$pval, scientific = TRUE), ")\n"
  )
}

report_text <- paste0(report_text,
  "\nTop 5 Cell Type Pairs by Mean Probability:\n"
)

for (i in 1:min(5, nrow(top_cell_pairs))) {
  pair <- top_cell_pairs[i, ]
  report_text <- paste0(report_text,
    "  ", i, ". ", pair$source, " → ", pair$target, 
    ": ", round(pair$mean_prob, 4), " (", pair$n_interactions, " interactions)\n"
  )
}

report_text <- paste0(report_text,
  "\nVisualization Files Created:\n",
  "  - Dot plots with different color palettes\n",
  "  - Bubble plot showing significance and probability\n",
  "  - Cell type communication heatmap\n"
)

if (exists("pathway_result") && !is.null(pathway_result)) {
  report_text <- paste0(report_text, "  - Pathway analysis plot\n")
}

report_text <- paste0(report_text,
  "  - Network summary plot\n\n",
  "Output Files:\n",
  "  - Communication matrix: communication_matrix_", JOB_ID, ".csv\n",
  "  - Cell type summary: celltype_communication_summary_", JOB_ID, ".csv\n",
  "  - Network statistics: network_statistics_", JOB_ID, ".csv\n"
)

if (exists("pathway_result") && !is.null(pathway_result)) {
  report_text <- paste0(report_text, "  - Pathway summary: pathway_summary_", JOB_ID, ".csv\n")
}

report_text <- paste0(report_text,
  "\nRecommendations:\n",
  "1. Focus on high-probability interactions for validation\n",
  "2. Investigate cell type pairs with strong communication\n",
  "3. Examine pathway-specific communication patterns\n",
  "4. Consider spatial context for communication validation\n",
  "5. Compare communication patterns across conditions\n\n",
  "Usage Examples:\n",
  "Rscript communication_visualization.R input_data.csv job_id\n",
  "Rscript communication_visualization.R cellchat_results.csv analysis_v1\n"
)

# Save report
report_file <- file.path(OUTPUT_DIR, paste0("visualization_report_", JOB_ID, ".txt"))
writeLines(report_text, report_file)

# Save session info
session_file <- file.path(OUTPUT_DIR, paste0("sessionInfo_", JOB_ID, ".txt"))
writeLines(capture.output(sessionInfo()), session_file)

# ============================================================================
# Final Summary
# ============================================================================

cat("\n", "="*60, "\n")
cat("Communication visualization analysis completed!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")
cat("\nSummary:\n")
cat("- Total interactions analyzed:", total_interactions, "\n")
cat("- Unique L-R pairs:", unique_lr_pairs, "\n")
cat("- Unique cell type pairs:", unique_cell_pairs, "\n")
cat("- Mean communication probability:", round(mean_prob, 4), "\n")

if (exists("top_interactions") && nrow(top_interactions) > 0) {
  cat("- Top interaction:", 
      paste(top_interactions$source[1], "→", top_interactions$target[1], 
            ":", top_interactions$lr_pair[1]), "\n")
}

cat("\nNext steps:\n")
cat("1. Review generated visualizations for biological insights\n")
cat("2. Validate key interactions with expression data\n")
cat("3. Compare communication patterns across samples\n")
cat("4. Investigate spatial distribution of communicating cells\n")

cat("\nAll visualizations and analysis results saved successfully!\n")