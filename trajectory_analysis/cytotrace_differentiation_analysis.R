#!/usr/bin/env Rscript

# ============================================================================
# CytoTRACE2 Differentiation Potential Analysis
# ============================================================================
# 
# This script performs cellular differentiation potential analysis using
# CytoTRACE2, which predicts the differentiation state of single cells
# based on transcriptional diversity and gene expression patterns.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE2)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(pheatmap)
  library(corrplot)
})

# Parameters
SPECIES <- "human"                    # "human" or "mouse"
SLOT_TYPE <- "counts"                 # "counts" or "data"
IS_SEURAT <- TRUE                     # Input is Seurat object
NCORES <- 4                          # Number of cores for parallel processing
BATCH_SIZE <- 25000                  # Batch size for large datasets
SMOOTH_BATCH_SIZE <- 1000            # Batch size for smoothing

# Plotting parameters
POINT_SIZE <- 0.5                    # Point size for plots
ALPHA <- 0.8                         # Transparency
COLOR_PALETTE <- "Spectral"          # Color palette for differentiation

# Input/Output paths
INPUT_SEURAT <- "../data/annotated_seurat_object.rds"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  INPUT_SEURAT <- args[1]
}
if (length(args) > 1) {
  SPECIES <- args[2]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to create differentiation color palette
create_differentiation_palette <- function(n_colors = 100) {
  # Create palette from more differentiated (blue) to less differentiated (red)
  if (COLOR_PALETTE == "Spectral") {
    colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(n_colors)
  } else if (COLOR_PALETTE == "RdYlBu") {
    colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(n_colors)
  } else {
    colors <- viridis(n_colors, option = "plasma")
  }
  return(colors)
}

# Function to calculate correlation with known markers
calculate_marker_correlations <- function(sobj, cytotrace_scores) {
  
  # Define differentiation markers for different systems
  differentiation_markers <- list(
    # Stemness/pluripotency markers (should be positively correlated with CytoTRACE2)
    stemness = c("NANOG", "POU5F1", "SOX2", "KLF4", "MYC", "LIN28A", "DPPA4", "DPPA2"),
    
    # Early differentiation markers (should be negatively correlated)
    early_diff = c("PAX6", "GATA1", "GATA2", "RUNX1", "HAND1", "HAND2", "TBX1", "TBX5"),
    
    # Mature/terminal differentiation markers (should be strongly negatively correlated)
    terminal_diff = c("ALB", "INS", "MYH6", "ACTC1", "NKX2-5", "TNNT2", "MBP", "SYP"),
    
    # Cell cycle markers (may correlate with differentiation potential)
    cell_cycle = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNA2")
  )
  
  correlations <- list()
  
  for (marker_type in names(differentiation_markers)) {
    markers <- differentiation_markers[[marker_type]]
    available_markers <- markers[markers %in% rownames(sobj)]
    
    if (length(available_markers) > 0) {
      # Calculate average expression of available markers
      if (length(available_markers) == 1) {
        marker_expression <- GetAssayData(sobj, slot = "data")[available_markers, ]
      } else {
        marker_expression <- colMeans(GetAssayData(sobj, slot = "data")[available_markers, ])
      }
      
      # Calculate correlation with CytoTRACE2 scores
      cor_result <- cor.test(cytotrace_scores, marker_expression, method = "spearman")
      
      correlations[[marker_type]] <- list(
        correlation = cor_result$estimate,
        p_value = cor_result$p.value,
        n_markers = length(available_markers),
        markers_used = available_markers
      )
    }
  }
  
  return(correlations)
}

# Function to analyze differentiation by cell type
analyze_differentiation_by_celltype <- function(sobj, cell_type_col = "celltype") {
  
  if (!cell_type_col %in% colnames(sobj@meta.data)) {
    cat("Warning: Cell type column", cell_type_col, "not found\\n")
    return(NULL)
  }
  
  # Calculate differentiation statistics by cell type
  diff_stats <- sobj@meta.data %>%
    group_by(!!sym(cell_type_col)) %>%
    summarise(
      n_cells = n(),
      mean_cytotrace = mean(CytoTRACE2_Relative, na.rm = TRUE),
      median_cytotrace = median(CytoTRACE2_Relative, na.rm = TRUE),
      sd_cytotrace = sd(CytoTRACE2_Relative, na.rm = TRUE),
      q25_cytotrace = quantile(CytoTRACE2_Relative, 0.25, na.rm = TRUE),
      q75_cytotrace = quantile(CytoTRACE2_Relative, 0.75, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_cytotrace))
  
  # Add differentiation categories
  diff_stats$differentiation_level <- cut(
    diff_stats$mean_cytotrace,
    breaks = c(0, 0.3, 0.7, 1.0),
    labels = c("Highly Differentiated", "Intermediate", "Poorly Differentiated"),
    include.lowest = TRUE
  )
  
  return(diff_stats)
}

# Function to identify trajectory-related genes
identify_differentiation_genes <- function(sobj, n_genes = 50) {
  
  # Find genes correlated with differentiation potential
  expr_data <- GetAssayData(sobj, slot = "data")
  cytotrace_scores <- sobj$CytoTRACE2_Relative
  
  # Calculate correlations for all genes
  gene_correlations <- apply(expr_data, 1, function(gene_expr) {
    if (sum(gene_expr > 0) < 10) return(NA)  # Skip lowly expressed genes
    
    cor_result <- cor.test(cytotrace_scores, gene_expr, method = "spearman")
    return(cor_result$estimate)
  })
  
  # Remove NA values
  gene_correlations <- gene_correlations[!is.na(gene_correlations)]
  
  # Get top positively and negatively correlated genes
  top_positive <- head(sort(gene_correlations, decreasing = TRUE), n_genes/2)
  top_negative <- head(sort(gene_correlations, decreasing = FALSE), n_genes/2)
  
  differentiation_genes <- list(
    stemness_associated = names(top_positive),
    differentiation_associated = names(top_negative),
    correlations = gene_correlations[c(names(top_positive), names(top_negative))]
  )
  
  return(differentiation_genes)
}

# ============================================================================
# Load and Prepare Data
# ============================================================================

cat("Loading Seurat object for CytoTRACE2 analysis...\\n")

if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded Seurat object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Check for UMAP coordinates
if (!"umap" %in% names(sobj@reductions)) {
  cat("UMAP not found. Computing UMAP...\\n")
  sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)
}

# Determine cell type column
cell_type_col <- "celltype"
if (!"celltype" %in% colnames(sobj@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "preliminary_celltype"
  } else if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "seurat_clusters"
  } else {
    cell_type_col <- NULL
    cat("Warning: No cell type annotation found\\n")
  }
}

if (!is.null(cell_type_col)) {
  cat("Using cell type column:", cell_type_col, "\\n")
}

# ============================================================================
# Run CytoTRACE2 Analysis
# ============================================================================

cat("Running CytoTRACE2 analysis...\\n")
cat("Species:", SPECIES, "\\n")
cat("Number of cores:", NCORES, "\\n")

# Check if CytoTRACE2 scores already exist
if ("CytoTRACE2_Relative" %in% colnames(sobj@meta.data)) {
  cat("CytoTRACE2 scores already exist. Using existing scores.\\n")
  cat("To recompute, remove the 'CytoTRACE2_Relative' column from metadata.\\n")
} else {
  
  # Prepare data for CytoTRACE2
  if (ncol(sobj) > BATCH_SIZE) {
    cat("Large dataset detected. Processing in batches...\\n")
    
    # Split data into batches
    n_batches <- ceiling(ncol(sobj) / BATCH_SIZE)
    batch_results <- list()
    
    for (i in 1:n_batches) {
      start_idx <- (i - 1) * BATCH_SIZE + 1
      end_idx <- min(i * BATCH_SIZE, ncol(sobj))
      
      cat("Processing batch", i, "/", n_batches, "(cells", start_idx, "-", end_idx, ")...\\n")
      
      # Extract batch
      sobj_batch <- sobj[, start_idx:end_idx]
      
      # Run CytoTRACE2 on batch
      tryCatch({
        cytotrace_result <- cytotrace2(
          sobj_batch,
          species = SPECIES,
          slot_type = SLOT_TYPE,
          is_seurat = IS_SEURAT,
          ncores = NCORES,
          batch_size = min(BATCH_SIZE, ncol(sobj_batch)),
          smooth_batch_size = min(SMOOTH_BATCH_SIZE, ncol(sobj_batch))
        )
        
        batch_results[[i]] <- cytotrace_result
        
      }, error = function(e) {
        cat("Error in batch", i, ":", e$message, "\\n")
        batch_results[[i]] <- NULL
      })
    }
    
    # Combine batch results
    valid_batches <- batch_results[!sapply(batch_results, is.null)]
    
    if (length(valid_batches) > 0) {
      # Combine CytoTRACE2 scores
      all_scores <- do.call(c, lapply(valid_batches, function(x) x$CytoTRACE2_Relative))
      sobj$CytoTRACE2_Relative <- all_scores[colnames(sobj)]
      
      cat("CytoTRACE2 analysis completed for", length(all_scores), "cells\\n")
    } else {
      stop("All batches failed. Please check data format and parameters.")
    }
    
  } else {
    # Run CytoTRACE2 on entire dataset
    cat("Running CytoTRACE2 on full dataset...\\n")
    
    tryCatch({
      cytotrace_result <- cytotrace2(
        sobj,
        species = SPECIES,
        slot_type = SLOT_TYPE,
        is_seurat = IS_SEURAT,
        ncores = NCORES,
        batch_size = BATCH_SIZE,
        smooth_batch_size = SMOOTH_BATCH_SIZE
      )
      
      # Add scores to Seurat object
      sobj$CytoTRACE2_Relative <- cytotrace_result$CytoTRACE2_Relative
      
      cat("CytoTRACE2 analysis completed successfully\\n")
      
    }, error = function(e) {
      stop("CytoTRACE2 analysis failed: ", e$message)
    })
  }
}

# Verify CytoTRACE2 scores
if (!"CytoTRACE2_Relative" %in% colnames(sobj@meta.data)) {
  stop("CytoTRACE2 scores not found in metadata")
}

cytotrace_range <- range(sobj$CytoTRACE2_Relative, na.rm = TRUE)
cat("CytoTRACE2 score range: [", round(cytotrace_range[1], 3), ", ", 
    round(cytotrace_range[2], 3), "]\\n")

# ============================================================================
# Visualization
# ============================================================================

cat("Creating differentiation potential visualizations...\\n")

# Create differentiation color palette
diff_colors <- create_differentiation_palette()

# 1. UMAP colored by CytoTRACE2 scores
umap_diff <- DimPlot(sobj, reduction = "umap", pt.size = POINT_SIZE) +
  geom_point(aes(color = sobj$CytoTRACE2_Relative), size = POINT_SIZE, alpha = ALPHA) +
  scale_color_gradientn(
    colors = diff_colors,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0.0\\n(More diff.)", "0.5", "1.0\\n(Less diff.)"),
    name = "Differentiation\\nPotential"
  ) +
  ggtitle("Cellular Differentiation Potential (CytoTRACE2)") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(file.path(PLOTS_DIR, "cytotrace2_umap_differentiation.pdf"), 
       umap_diff, width = 10, height = 8)

# 2. UMAP colored by cell type (if available)
if (!is.null(cell_type_col)) {
  umap_celltype <- DimPlot(sobj, reduction = "umap", group.by = cell_type_col, 
                          pt.size = POINT_SIZE, label = TRUE, repel = TRUE) +
    ggtitle("Cell Types") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    coord_fixed()
  
  # Combine plots
  combined_umap <- umap_celltype | umap_diff
  ggsave(file.path(PLOTS_DIR, "cytotrace2_celltype_comparison.pdf"), 
         combined_umap, width = 20, height = 8)
}

# 3. Differentiation distribution by cell type
if (!is.null(cell_type_col)) {
  diff_violin <- ggplot(sobj@meta.data, aes_string(x = cell_type_col, y = "CytoTRACE2_Relative")) +
    geom_violin(aes_string(fill = cell_type_col), alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    scale_fill_discrete(guide = "none") +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0.0\n(More diff.)", "0.5", "1.0\n(Less diff.)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    ) +
    labs(
      title = "Differentiation Potential Distribution by Cell Type",
      y = "CytoTRACE2 Score"
    )
  
  ggsave(file.path(PLOTS_DIR, "cytotrace2_distribution_by_celltype.pdf"), 
         diff_violin, width = 12, height = 6)
}

# 4. Differentiation histogram
diff_hist <- ggplot(sobj@meta.data, aes(x = CytoTRACE2_Relative)) +
  geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "black") +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0.0\n(More diff.)", "0.5", "1.0\n(Less diff.)")
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of Differentiation Potential",
    x = "CytoTRACE2 Score",
    y = "Number of Cells"
  )

ggsave(file.path(PLOTS_DIR, "cytotrace2_score_distribution.pdf"), 
       diff_hist, width = 8, height = 6)

# ============================================================================
# Statistical Analysis
# ============================================================================

cat("Performing statistical analysis...\\n")

# Analyze differentiation by cell type
if (!is.null(cell_type_col)) {
  diff_by_celltype <- analyze_differentiation_by_celltype(sobj, cell_type_col)
  
  if (!is.null(diff_by_celltype)) {
    write.csv(diff_by_celltype, file.path(OUTPUT_DIR, "differentiation_by_celltype.csv"), 
              row.names = FALSE)
    
    cat("Cell types ranked by differentiation potential (least to most differentiated):\\n")
    print(diff_by_celltype[, c(cell_type_col, "mean_cytotrace", "differentiation_level")])
  }
}

# Calculate correlations with known markers
marker_correlations <- calculate_marker_correlations(sobj, sobj$CytoTRACE2_Relative)

if (length(marker_correlations) > 0) {
  # Create correlation summary
  cor_summary <- data.frame(
    marker_type = names(marker_correlations),
    correlation = sapply(marker_correlations, function(x) x$correlation),
    p_value = sapply(marker_correlations, function(x) x$p_value),
    n_markers = sapply(marker_correlations, function(x) x$n_markers),
    stringsAsFactors = FALSE
  )
  
  write.csv(cor_summary, file.path(OUTPUT_DIR, "marker_correlations.csv"), 
            row.names = FALSE)
  
  cat("\\nCorrelations with known marker sets:\\n")
  print(cor_summary)
}

# Identify differentiation-associated genes
diff_genes <- identify_differentiation_genes(sobj, n_genes = 100)

if (length(diff_genes$stemness_associated) > 0) {
  # Save gene lists
  write.table(diff_genes$stemness_associated, 
              file.path(OUTPUT_DIR, "stemness_associated_genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  write.table(diff_genes$differentiation_associated, 
              file.path(OUTPUT_DIR, "differentiation_associated_genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Save correlations
  gene_cors <- data.frame(
    gene = names(diff_genes$correlations),
    correlation = diff_genes$correlations,
    type = ifelse(diff_genes$correlations > 0, "Stemness", "Differentiation"),
    stringsAsFactors = FALSE
  )
  
  write.csv(gene_cors, file.path(OUTPUT_DIR, "differentiation_gene_correlations.csv"), 
            row.names = FALSE)
  
  cat("\\nTop stemness-associated genes:", 
      paste(head(diff_genes$stemness_associated, 5), collapse = ", "), "\\n")
  cat("Top differentiation-associated genes:", 
      paste(head(diff_genes$differentiation_associated, 5), collapse = ", "), "\\n")
}

# ============================================================================
# Advanced Visualizations
# ============================================================================

cat("Creating advanced visualizations...\\n")

# 5. Feature plots for top differentiation genes
if (length(diff_genes$stemness_associated) >= 4) {
  top_stemness_genes <- head(diff_genes$stemness_associated, 4)
  
  stemness_plots <- list()
  for (i in seq_along(top_stemness_genes)) {
    stemness_plots[[i]] <- FeaturePlot(sobj, features = top_stemness_genes[i], 
                                      pt.size = POINT_SIZE/2, order = TRUE) +
      theme(axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), plot.title = element_text(size = 10))
  }
  
  combined_stemness <- wrap_plots(stemness_plots, ncol = 2)
  ggsave(file.path(PLOTS_DIR, "top_stemness_genes_features.pdf"), 
         combined_stemness, width = 12, height = 10)
}

if (length(diff_genes$differentiation_associated) >= 4) {
  top_diff_genes <- head(diff_genes$differentiation_associated, 4)
  
  diff_plots <- list()
  for (i in seq_along(top_diff_genes)) {
    diff_plots[[i]] <- FeaturePlot(sobj, features = top_diff_genes[i], 
                                  pt.size = POINT_SIZE/2, order = TRUE) +
      theme(axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), plot.title = element_text(size = 10))
  }
  
  combined_diff <- wrap_plots(diff_plots, ncol = 2)
  ggsave(file.path(PLOTS_DIR, "top_differentiation_genes_features.pdf"), 
         combined_diff, width = 12, height = 10)
}

# 6. Correlation heatmap between CytoTRACE2 and selected genes
if (length(diff_genes$correlations) > 0) {
  # Select top genes for heatmap
  top_cor_genes <- names(sort(abs(diff_genes$correlations), decreasing = TRUE))[1:min(20, length(diff_genes$correlations))]
  
  if (length(top_cor_genes) > 1) {
    # Get expression data for selected genes
    expr_data <- GetAssayData(sobj, slot = "data")[top_cor_genes, ]
    
    # Calculate correlation matrix
    cor_matrix <- cor(t(rbind(sobj$CytoTRACE2_Relative, expr_data)), method = "spearman")
    cor_matrix <- cor_matrix[-1, 1, drop = FALSE]  # Remove self-correlation
    
    # Create heatmap
    pdf(file.path(PLOTS_DIR, "cytotrace2_gene_correlations.pdf"), width = 8, height = 10)
    pheatmap(
      cor_matrix,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      main = "Gene Correlations with CytoTRACE2",
      fontsize = 10
    )
    dev.off()
  }
}

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating comprehensive report...\\n")

# Calculate summary statistics
n_cells <- ncol(sobj)
cytotrace_mean <- mean(sobj$CytoTRACE2_Relative, na.rm = TRUE)
cytotrace_sd <- sd(sobj$CytoTRACE2_Relative, na.rm = TRUE)

# Create report
report_text <- paste0(
  "CytoTRACE2 Differentiation Analysis Report\\n",
  "=========================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Species: ", SPECIES, "\\n\\n",
  "Data Summary:\\n",
  "  - Total cells: ", n_cells, "\\n",
  "  - CytoTRACE2 score range: [", round(cytotrace_range[1], 3), ", ", 
    round(cytotrace_range[2], 3), "]\\n",
  "  - Mean CytoTRACE2 score: ", round(cytotrace_mean, 3), "\\n",
  "  - Standard deviation: ", round(cytotrace_sd, 3), "\\n\\n",
  "Analysis Parameters:\\n",
  "  - Species: ", SPECIES, "\\n",
  "  - Slot type: ", SLOT_TYPE, "\\n",
  "  - Number of cores: ", NCORES, "\\n",
  "  - Batch size: ", BATCH_SIZE, "\\n\\n"
)

# Add cell type analysis if available
if (!is.null(cell_type_col) && !is.null(diff_by_celltype)) {
  most_stem <- diff_by_celltype[[cell_type_col]][1]
  most_diff <- diff_by_celltype[[cell_type_col]][nrow(diff_by_celltype)]
  
  report_text <- paste0(report_text,
    "Cell Type Analysis:\\n",
    "  - Most stem-like: ", most_stem, " (", 
         round(diff_by_celltype$mean_cytotrace[1], 3), ")\\n",
    "  - Most differentiated: ", most_diff, " (", 
         round(diff_by_celltype$mean_cytotrace[nrow(diff_by_celltype)], 3), ")\\n",
    "  - Cell types analyzed: ", nrow(diff_by_celltype), "\\n\\n"
  )
}

# Add marker correlation results
if (length(marker_correlations) > 0) {
  report_text <- paste0(report_text,
    "Marker Validation:\\n"
  )
  
  for (marker_type in names(marker_correlations)) {
    cor_val <- marker_correlations[[marker_type]]$correlation
    p_val <- marker_correlations[[marker_type]]$p_value
    
    report_text <- paste0(report_text,
      "  - ", marker_type, ": r = ", round(cor_val, 3), 
      " (p = ", format(p_val, scientific = TRUE, digits = 2), ")\\n"
    )
  }
  report_text <- paste0(report_text, "\\n")
}

# Add gene analysis results
if (length(diff_genes$stemness_associated) > 0) {
  report_text <- paste0(report_text,
    "Differentiation-Associated Genes:\\n",
    "  - Stemness-associated genes: ", length(diff_genes$stemness_associated), "\\n",
    "  - Differentiation-associated genes: ", length(diff_genes$differentiation_associated), "\\n\\n"
  )
}

report_text <- paste0(report_text,
  "Output Files:\\n",
  "  - Processed Seurat object: cytotrace2_analysis.rds\\n"
)

if (!is.null(diff_by_celltype)) {
  report_text <- paste0(report_text,
    "  - Differentiation by cell type: differentiation_by_celltype.csv\\n"
  )
}

if (length(marker_correlations) > 0) {
  report_text <- paste0(report_text,
    "  - Marker correlations: marker_correlations.csv\\n"
  )
}

if (length(diff_genes$correlations) > 0) {
  report_text <- paste0(report_text,
    "  - Gene correlations: differentiation_gene_correlations.csv\\n",
    "  - Stemness genes: stemness_associated_genes.txt\\n",
    "  - Differentiation genes: differentiation_associated_genes.txt\\n"
  )
}

report_text <- paste0(report_text,
  "\\nVisualization Files:\\n",
  "  - UMAP differentiation: cytotrace2_umap_differentiation.pdf\\n",
  "  - Score distribution: cytotrace2_score_distribution.pdf\\n"
)

if (!is.null(cell_type_col)) {
  report_text <- paste0(report_text,
    "  - Cell type comparison: cytotrace2_celltype_comparison.pdf\\n",
    "  - Distribution by cell type: cytotrace2_distribution_by_celltype.pdf\\n"
  )
}

report_text <- paste0(report_text,
  "\\nRecommendations:\\n",
  "1. Validate stem-like populations with known stemness markers\\n",
  "2. Investigate highly differentiated cells for terminal markers\\n",
  "3. Use CytoTRACE2 scores for trajectory analysis root selection\\n",
  "4. Compare with other differentiation potential methods\\n",
  "5. Analyze differentiation-associated genes for regulatory mechanisms\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "cytotrace2_analysis_report.txt"))

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\\n")

# Save processed Seurat object
saveRDS(sobj, file.path(OUTPUT_DIR, "cytotrace2_analysis.rds"))

# Save CytoTRACE2 scores
cytotrace_scores <- data.frame(
  cell_id = colnames(sobj),
  CytoTRACE2_Relative = sobj$CytoTRACE2_Relative,
  stringsAsFactors = FALSE
)

if (!is.null(cell_type_col)) {
  cytotrace_scores[[cell_type_col]] <- sobj@meta.data[[cell_type_col]]
}

write.csv(cytotrace_scores, file.path(OUTPUT_DIR, "cytotrace2_scores.csv"), 
          row.names = FALSE)

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_cytotrace2_analysis.txt"))

# ============================================================================
# Final Summary
# ============================================================================

cat("\\nCytoTRACE2 differentiation analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nSummary:\\n")
cat("- Cells analyzed:", n_cells, "\\n")
cat("- CytoTRACE2 range: [", round(cytotrace_range[1], 3), ", ", 
    round(cytotrace_range[2], 3), "]\\n")
cat("- Mean differentiation potential:", round(cytotrace_mean, 3), "\\n")

if (!is.null(diff_by_celltype)) {
  cat("- Cell types analyzed:", nrow(diff_by_celltype), "\\n")
}

if (length(marker_correlations) > 0) {
  cat("- Marker correlations calculated:", length(marker_correlations), "\\n")
}

if (length(diff_genes$correlations) > 0) {
  cat("- Differentiation genes identified:", length(diff_genes$correlations), "\\n")
}

cat("\\nNext steps:\\n")
cat("1. Validate CytoTRACE2 predictions with experimental data\\n")
cat("2. Use scores for trajectory analysis root identification\\n")
cat("3. Investigate cell type-specific differentiation patterns\\n")
cat("4. Analyze differentiation-associated gene networks\\n")

cat("\\nUsage examples:\\n")
cat("Rscript cytotrace_differentiation_analysis.R <input_seurat.rds> <species>\\n")
cat("Rscript cytotrace_differentiation_analysis.R data.rds human\\n")