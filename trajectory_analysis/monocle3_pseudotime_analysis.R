#!/usr/bin/env Rscript

# ============================================================================
# Monocle3 Pseudotime and Trajectory Analysis
# ============================================================================
# 
# This script performs pseudotime analysis using Monocle3 to infer developmental
# trajectories and identify branching points in cell differentiation. It 
# integrates with CytoTRACE2 for differentiation potential analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(RColorBrewer)
  library(viridis)
  library(scales)
})

# Parameters
CLUSTERING_RESOLUTION <- 1e-3          # Low resolution for trajectory clustering
MINIMAL_BRANCH_LEN_RANGE <- c(10, 50)  # Range for branch length optimization
NCENTER_RANGE <- c(50, 1000)           # Range for ncenter optimization
N_CORES <- 4                           # Number of cores for parallel processing
USE_PARTITION <- TRUE                  # Whether to use partition-based trajectory learning
CLOSE_LOOP <- FALSE                    # Whether to allow closed loops in trajectories

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

# ============================================================================
# Helper Functions
# ============================================================================

# Function to create differentiation color palette
create_differentiation_palette <- function() {
  # Spectral palette from more differentiated (blue) to less differentiated (red)
  colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
  return(colors)
}

# Function to optimize trajectory parameters
optimize_trajectory_parameters <- function(cds, output_prefix, test_minimal_branch = TRUE, test_ncenter = TRUE) {
  
  param_results <- data.frame(
    parameter = character(),
    value = numeric(),
    n_branches = numeric(),
    n_leaves = numeric(),
    connectivity = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Test minimal_branch_len if requested
  if (test_minimal_branch) {
    cat("Optimizing minimal_branch_len parameter...\\n")
    
    for (len in seq(MINIMAL_BRANCH_LEN_RANGE[1], MINIMAL_BRANCH_LEN_RANGE[2], by = 5)) {
      tryCatch({
        cds_temp <- learn_graph(cds, 
                               use_partition = USE_PARTITION, 
                               verbose = FALSE, 
                               close_loop = CLOSE_LOOP,
                               learn_graph_control = list(minimal_branch_len = len))
        
        # Calculate trajectory statistics
        n_branches <- length(unique(principal_graph(cds_temp)$UMAP$branch_id))
        n_leaves <- sum(principal_graph(cds_temp)$UMAP$is_leaf)
        
        # Estimate connectivity (simplified)
        connectivity <- length(principal_graph(cds_temp)$UMAP$edge_links) / 
                       nrow(principal_graph(cds_temp)$UMAP)
        
        param_results <- rbind(param_results, data.frame(
          parameter = "minimal_branch_len",
          value = len,
          n_branches = n_branches,
          n_leaves = n_leaves,
          connectivity = connectivity,
          stringsAsFactors = FALSE
        ))
        
        # Save diagnostic plot
        p <- plot_cells(cds_temp,
                       color_cells_by = "seurat_clusters",
                       label_groups_by_cluster = FALSE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE) +
          ggtitle(paste("Branch Length =", len))
        
        ggsave(file.path(PLOTS_DIR, paste0(output_prefix, "_branch_len_", len, ".pdf")), 
               p, width = 8, height = 6)
        
      }, error = function(e) {
        cat("Error with minimal_branch_len =", len, ":", e$message, "\\n")
      })
    }
  }
  
  # Test ncenter if requested
  if (test_ncenter) {
    cat("Optimizing ncenter parameter...\\n")
    
    for (ncenter_val in seq(NCENTER_RANGE[1], NCENTER_RANGE[2], by = 50)) {
      tryCatch({
        cds_temp <- learn_graph(cds, 
                               use_partition = USE_PARTITION, 
                               verbose = FALSE, 
                               close_loop = CLOSE_LOOP,
                               learn_graph_control = list(ncenter = ncenter_val))
        
        # Calculate trajectory statistics
        n_branches <- length(unique(principal_graph(cds_temp)$UMAP$branch_id))
        n_leaves <- sum(principal_graph(cds_temp)$UMAP$is_leaf)
        
        # Estimate connectivity
        connectivity <- length(principal_graph(cds_temp)$UMAP$edge_links) / 
                       nrow(principal_graph(cds_temp)$UMAP)
        
        param_results <- rbind(param_results, data.frame(
          parameter = "ncenter",
          value = ncenter_val,
          n_branches = n_branches,
          n_leaves = n_leaves,
          connectivity = connectivity,
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        cat("Error with ncenter =", ncenter_val, ":", e$message, "\\n")
      })
    }
  }
  
  return(param_results)
}

# Function to analyze trajectory branches
analyze_trajectory_branches <- function(cds) {
  
  # Get principal graph information
  pg <- principal_graph(cds)$UMAP
  
  # Identify branch points and leaves
  branch_points <- pg[pg$branch_points, ]
  leaves <- pg[pg$is_leaf, ]
  
  # Calculate branch statistics
  branch_stats <- list(
    n_branch_points = nrow(branch_points),
    n_leaves = nrow(leaves),
    n_edges = nrow(pg),
    avg_branch_length = mean(pg$branch_length, na.rm = TRUE)
  )
  
  return(branch_stats)
}

# Function to identify root cells automatically
identify_root_cells <- function(cds, method = "earliest_pseudotime") {
  
  if (method == "earliest_pseudotime") {
    # Use cells with highest differentiation potential (if available)
    if ("CytoTRACE2_Relative" %in% colnames(colData(cds))) {
      # Cells with highest CytoTRACE2 scores (least differentiated)
      cytotrace_scores <- colData(cds)$CytoTRACE2_Relative
      root_cells <- colnames(cds)[which(cytotrace_scores > quantile(cytotrace_scores, 0.9, na.rm = TRUE))]
    } else {
      # Fallback: use cluster with highest gene expression diversity
      cluster_diversity <- aggregate(rowSums(exprs(cds) > 0), 
                                   by = list(clusters(cds)), 
                                   FUN = mean)
      root_cluster <- cluster_diversity$Group.1[which.max(cluster_diversity$x)]
      root_cells <- colnames(cds)[clusters(cds) == root_cluster]
    }
  } else if (method == "centroid") {
    # Use cells closest to centroid in embedding space
    embedding <- reducedDims(cds)$UMAP
    centroid <- colMeans(embedding)
    distances <- rowSums((embedding - matrix(centroid, nrow = nrow(embedding), 
                                           ncol = ncol(embedding), byrow = TRUE))^2)
    root_cells <- colnames(cds)[which(distances == min(distances))]
  }
  
  return(root_cells)
}

# ============================================================================
# Load and Prepare Data
# ============================================================================

cat("Loading Seurat object...\\n")

if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded Seurat object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Check for required reductions and metadata
if (!"umap" %in% names(sobj@reductions)) {
  cat("UMAP not found. Computing UMAP...\\n")
  sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)
}

# Check for cell type annotations
cell_type_col <- "celltype"
if (!"celltype" %in% colnames(sobj@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "preliminary_celltype"
  } else if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "seurat_clusters"
    sobj$celltype <- sobj$seurat_clusters
  } else {
    stop("No suitable cell type annotation found")
  }
}

cat("Using cell type column:", cell_type_col, "\\n")

# ============================================================================
# Convert to Monocle3 Format
# ============================================================================

cat("Converting Seurat object to Monocle3 cell_data_set...\\n")

# Convert to cell_data_set
cds <- as.cell_data_set(sobj)

# Transfer additional metadata if available
if ("CytoTRACE2_Relative" %in% colnames(sobj@meta.data)) {
  colData(cds)$CytoTRACE2_Relative <- sobj$CytoTRACE2_Relative
  cat("CytoTRACE2 differentiation scores transferred\\n")
}

# Transfer UMAP coordinates
if ("umap" %in% names(sobj@reductions)) {
  reducedDims(cds)$UMAP <- Embeddings(sobj, "umap")
  cat("UMAP coordinates transferred\\n")
}

cat("Cell_data_set created with", ncol(cds), "cells\\n")

# ============================================================================
# Clustering for Trajectory Analysis
# ============================================================================

cat("Performing clustering for trajectory analysis...\\n")

# Cluster cells with low resolution for trajectory learning
cds <- cluster_cells(cds, resolution = CLUSTERING_RESOLUTION, verbose = FALSE)

# Create cluster and partition plots
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE) +
  ggtitle("Monocle3 Clusters")

p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) +
  ggtitle("Monocle3 Partitions")

# Save clustering results
cluster_plots <- p1 | p2
ggsave(file.path(PLOTS_DIR, "monocle3_clustering_overview.pdf"), 
       cluster_plots, width = 14, height = 6)

cat("Identified", length(unique(clusters(cds))), "clusters and", 
    length(unique(partitions(cds))), "partitions\\n")

# ============================================================================
# Trajectory Learning and Optimization
# ============================================================================

cat("Learning trajectory graph...\\n")

# Option 1: Parameter optimization (can be time-consuming)
OPTIMIZE_PARAMETERS <- FALSE  # Set to TRUE for parameter optimization

if (OPTIMIZE_PARAMETERS) {
  cat("Optimizing trajectory parameters...\\n")
  param_results <- optimize_trajectory_parameters(cds, "trajectory_optimization", 
                                                  test_minimal_branch = TRUE, 
                                                  test_ncenter = FALSE)
  
  # Save parameter optimization results
  write.csv(param_results, file.path(OUTPUT_DIR, "trajectory_parameter_optimization.csv"), 
            row.names = FALSE)
  
  # Select optimal parameters based on connectivity and branch structure
  optimal_params <- param_results %>%
    filter(parameter == "minimal_branch_len") %>%
    filter(connectivity > quantile(connectivity, 0.3) & connectivity < quantile(connectivity, 0.7)) %>%
    filter(n_branches >= 2 & n_branches <= 10) %>%
    arrange(desc(connectivity)) %>%
    slice(1)
  
  if (nrow(optimal_params) > 0) {
    optimal_branch_len <- optimal_params$value
    cat("Using optimized minimal_branch_len:", optimal_branch_len, "\\n")
    
    cds <- learn_graph(cds, 
                      use_partition = USE_PARTITION, 
                      verbose = FALSE, 
                      close_loop = CLOSE_LOOP,
                      learn_graph_control = list(minimal_branch_len = optimal_branch_len))
  } else {
    cat("Using default parameters for trajectory learning\\n")
    cds <- learn_graph(cds, use_partition = USE_PARTITION, verbose = FALSE, close_loop = CLOSE_LOOP)
  }
} else {
  # Option 2: Use default parameters
  cat("Using default parameters for trajectory learning\\n")
  cds <- learn_graph(cds, use_partition = USE_PARTITION, verbose = FALSE, close_loop = CLOSE_LOOP)
}

# Analyze trajectory structure
branch_stats <- analyze_trajectory_branches(cds)
cat("Trajectory statistics:\\n")
cat("  - Branch points:", branch_stats$n_branch_points, "\\n")
cat("  - Leaves:", branch_stats$n_leaves, "\\n")
cat("  - Edges:", branch_stats$n_edges, "\\n")
cat("  - Average branch length:", round(branch_stats$avg_branch_length, 2), "\\n")

# ============================================================================
# Visualize Trajectory
# ============================================================================

cat("Creating trajectory visualizations...\\n")

# Create differentiation color palette
diff_colors <- create_differentiation_palette()

# 1. Basic trajectory plot by cell type
trajectory_celltype <- plot_cells(cds,
                                 color_cells_by = cell_type_col,
                                 label_groups_by_cluster = FALSE,
                                 label_leaves = FALSE,
                                 label_branch_points = FALSE,
                                 label_roots = FALSE) +
  ggtitle("Developmental Trajectory by Cell Type") +
  theme(legend.position = "bottom")

ggsave(file.path(PLOTS_DIR, "trajectory_by_celltype.pdf"), 
       trajectory_celltype, width = 10, height = 8)

# 2. Trajectory with CytoTRACE2 differentiation scores (if available)
if ("CytoTRACE2_Relative" %in% colnames(colData(cds))) {
  trajectory_diff <- plot_cells(cds,
                               color_cells_by = "CytoTRACE2_Relative",
                               label_groups_by_cluster = FALSE,
                               label_leaves = FALSE,
                               label_branch_points = FALSE,
                               label_roots = FALSE) +
    scale_color_gradientn(
      colors = diff_colors,
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("0.0 (More diff.)", "1.0 (Less diff.)"),
      name = "Differentiation\\nPotential"
    ) +
    ggtitle("Trajectory Colored by Differentiation Potential") +
    coord_fixed()
  
  ggsave(file.path(PLOTS_DIR, "trajectory_by_differentiation.pdf"), 
         trajectory_diff, width = 10, height = 8)
}

# 3. Trajectory with clusters
trajectory_clusters <- plot_cells(cds,
                                 color_cells_by = "cluster",
                                 label_groups_by_cluster = TRUE,
                                 label_leaves = TRUE,
                                 label_branch_points = TRUE,
                                 label_roots = FALSE) +
  ggtitle("Trajectory with Cluster Labels")

ggsave(file.path(PLOTS_DIR, "trajectory_with_labels.pdf"), 
       trajectory_clusters, width = 10, height = 8)

# ============================================================================
# Pseudotime Analysis
# ============================================================================

cat("Performing pseudotime analysis...\\n")

# Identify root cells automatically
if ("CytoTRACE2_Relative" %in% colnames(colData(cds))) {
  root_cells <- identify_root_cells(cds, method = "earliest_pseudotime")
  cat("Identified", length(root_cells), "root cells based on differentiation potential\\n")
} else {
  root_cells <- identify_root_cells(cds, method = "centroid")
  cat("Identified", length(root_cells), "root cells based on centroid method\\n")
}

# Order cells in pseudotime
if (length(root_cells) > 0) {
  cds <- order_cells(cds, root_cells = root_cells)
  
  # Check if pseudotime was calculated successfully
  if (all(is.finite(pseudotime(cds)))) {
    cat("Pseudotime calculated successfully\\n")
    
    # Pseudotime visualization
    pseudotime_plot <- plot_cells(cds,
                                 color_cells_by = "pseudotime",
                                 label_groups_by_cluster = FALSE,
                                 label_leaves = FALSE,
                                 label_branch_points = FALSE,
                                 label_roots = TRUE,
                                 trajectory_graph_color = "grey60") +
      scale_color_viridis_c(name = "Pseudotime") +
      ggtitle("Developmental Pseudotime")
    
    ggsave(file.path(PLOTS_DIR, "trajectory_pseudotime.pdf"), 
           pseudotime_plot, width = 10, height = 8)
    
    # Pseudotime by cell type
    pseudotime_by_celltype <- ggplot(as.data.frame(colData(cds)), 
                                    aes_string(x = cell_type_col, y = "pseudotime")) +
      geom_violin(aes_string(fill = cell_type_col), alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Pseudotime Distribution by Cell Type",
           x = "Cell Type", y = "Pseudotime")
    
    ggsave(file.path(PLOTS_DIR, "pseudotime_by_celltype.pdf"), 
           pseudotime_by_celltype, width = 12, height = 6)
    
    # Save pseudotime values
    pseudotime_df <- data.frame(
      cell_id = colnames(cds),
      pseudotime = pseudotime(cds),
      cell_type = colData(cds)[[cell_type_col]],
      cluster = clusters(cds),
      partition = partitions(cds),
      stringsAsFactors = FALSE
    )
    
    write.csv(pseudotime_df, file.path(OUTPUT_DIR, "pseudotime_values.csv"), 
              row.names = FALSE)
    
  } else {
    cat("Warning: Pseudotime calculation failed. Check trajectory structure.\\n")
  }
} else {
  cat("Warning: No root cells identified. Skipping pseudotime analysis.\\n")
}

# ============================================================================
# Gene Expression Along Pseudotime
# ============================================================================

if (all(is.finite(pseudotime(cds)))) {
  cat("Analyzing gene expression along pseudotime...\\n")
  
  # Find genes that change along pseudotime
  pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = N_CORES)
  
  # Filter significant genes
  pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
  
  if (length(pr_deg_ids) > 0) {
    cat("Found", length(pr_deg_ids), "genes significantly changing along pseudotime\\n")
    
    # Save pseudotime DEG results
    write.csv(pr_test_res, file.path(OUTPUT_DIR, "pseudotime_deg_results.csv"), 
              row.names = TRUE)
    
    # Plot top changing genes
    if (length(pr_deg_ids) >= 9) {
      top_genes <- head(pr_deg_ids[order(pr_test_res[pr_deg_ids, "morans_test_statistic"], 
                                        decreasing = TRUE)], 9)
      
      gene_plots <- plot_genes_in_pseudotime(cds[top_genes,],
                                            color_cells_by = cell_type_col,
                                            min_expr = 0.5,
                                            ncol = 3)
      
      ggsave(file.path(PLOTS_DIR, "top_genes_pseudotime.pdf"), 
             gene_plots, width = 15, height = 12)
    }
  } else {
    cat("No genes found significantly changing along pseudotime\\n")
  }
}

# ============================================================================
# Branch Analysis
# ============================================================================

if (branch_stats$n_branch_points > 0) {
  cat("Analyzing trajectory branches...\\n")
  
  # Branch-specific gene expression (requires manual selection of branches)
  # This would typically involve choose_graph_segments() interactively
  
  # For automated analysis, we can analyze cell fate decisions
  if (all(is.finite(pseudotime(cds)))) {
    
    # Identify cells at branch points (high pseudotime variance in neighborhood)
    # This is a simplified approach - more sophisticated methods available
    
    branch_cells_analysis <- data.frame(
      cell_id = colnames(cds),
      pseudotime = pseudotime(cds),
      cell_type = colData(cds)[[cell_type_col]],
      is_branch_cell = FALSE,
      stringsAsFactors = FALSE
    )
    
    # Mark cells near branch points based on trajectory graph
    # (This is a simplified implementation)
    pseudotime_range <- range(pseudotime(cds), na.rm = TRUE)
    branch_regions <- seq(pseudotime_range[1], pseudotime_range[2], length.out = 5)
    
    for (i in 2:(length(branch_regions)-1)) {
      region_cells <- which(abs(pseudotime(cds) - branch_regions[i]) < 
                           diff(pseudotime_range) * 0.1)
      if (length(region_cells) > 10) {
        region_celltypes <- table(colData(cds)[[cell_type_col]][region_cells])
        if (length(region_celltypes) >= 2) {
          branch_cells_analysis$is_branch_cell[region_cells] <- TRUE
        }
      }
    }
    
    # Save branch analysis
    write.csv(branch_cells_analysis, file.path(OUTPUT_DIR, "branch_cells_analysis.csv"), 
              row.names = FALSE)
    
    n_branch_cells <- sum(branch_cells_analysis$is_branch_cell)
    cat("Identified", n_branch_cells, "potential branch cells\\n")
  }
}

# ============================================================================
# Save Results and Generate Report
# ============================================================================

cat("Saving results and generating report...\\n")

# Save processed cell_data_set
saveRDS(cds, file.path(OUTPUT_DIR, "monocle3_trajectory_analysis.rds"))

# Generate comprehensive report
n_cells <- ncol(cds)
n_clusters <- length(unique(clusters(cds)))
n_partitions <- length(unique(partitions(cds)))

report_text <- paste0(
  "Monocle3 Trajectory Analysis Report\\n",
  "==================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n\\n",
  "Data Summary:\\n",
  "  - Total cells: ", n_cells, "\\n",
  "  - Monocle3 clusters: ", n_clusters, "\\n",
  "  - Partitions: ", n_partitions, "\\n",
  "  - Cell types: ", length(unique(colData(cds)[[cell_type_col]])), "\\n\\n",
  "Trajectory Structure:\\n",
  "  - Branch points: ", branch_stats$n_branch_points, "\\n",
  "  - Leaves: ", branch_stats$n_leaves, "\\n",
  "  - Edges: ", branch_stats$n_edges, "\\n",
  "  - Average branch length: ", round(branch_stats$avg_branch_length, 2), "\\n\\n",
  "Analysis Parameters:\\n",
  "  - Clustering resolution: ", CLUSTERING_RESOLUTION, "\\n",
  "  - Use partition: ", USE_PARTITION, "\\n",
  "  - Close loops: ", CLOSE_LOOP, "\\n",
  "  - Parameter optimization: ", OPTIMIZE_PARAMETERS, "\\n\\n"
)

if (all(is.finite(pseudotime(cds)))) {
  pseudotime_range <- range(pseudotime(cds), na.rm = TRUE)
  report_text <- paste0(report_text,
    "Pseudotime Analysis:\\n",
    "  - Pseudotime calculated: Yes\\n",
    "  - Root cells identified: ", length(root_cells), "\\n",
    "  - Pseudotime range: [", round(pseudotime_range[1], 2), ", ", 
         round(pseudotime_range[2], 2), "]\\n"
  )
  
  if (exists("pr_deg_ids")) {
    report_text <- paste0(report_text,
      "  - Pseudotime DEGs: ", length(pr_deg_ids), "\\n"
    )
  }
} else {
  report_text <- paste0(report_text,
    "Pseudotime Analysis:\\n",
    "  - Pseudotime calculated: No\\n",
    "  - Issue: Failed to identify suitable trajectory structure\\n"
  )
}

report_text <- paste0(report_text,
  "\\nOutput Files:\\n",
  "  - Processed CDS: monocle3_trajectory_analysis.rds\\n",
  "  - Pseudotime values: pseudotime_values.csv\\n"
)

if (exists("pr_test_res")) {
  report_text <- paste0(report_text,
    "  - Pseudotime DEGs: pseudotime_deg_results.csv\\n"
  )
}

if (exists("branch_cells_analysis")) {
  report_text <- paste0(report_text,
    "  - Branch analysis: branch_cells_analysis.csv\\n"
  )
}

report_text <- paste0(report_text,
  "\\nVisualization Files:\\n",
  "  - trajectory_by_celltype.pdf: Main trajectory plot\\n",
  "  - trajectory_pseudotime.pdf: Pseudotime visualization\\n",
  "  - pseudotime_by_celltype.pdf: Pseudotime distribution\\n",
  "  - monocle3_clustering_overview.pdf: Clustering results\\n\\n",
  "Recommendations:\\n",
  "1. Validate trajectory structure with biological knowledge\\n",
  "2. Investigate branch points for cell fate decisions\\n",
  "3. Analyze pseudotime DEGs for developmental programs\\n",
  "4. Consider trajectory-based subclustering for refinement\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "monocle3_trajectory_report.txt"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_monocle3_trajectory.txt"))

# ============================================================================
# Final Summary
# ============================================================================

cat("\\nMonocle3 trajectory analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nTrajectory summary:\\n")
cat("- Cells analyzed:", n_cells, "\\n")
cat("- Clusters identified:", n_clusters, "\\n")
cat("- Branch points:", branch_stats$n_branch_points, "\\n")

if (all(is.finite(pseudotime(cds)))) {
  cat("- Pseudotime range: [", round(pseudotime_range[1], 2), ", ", 
      round(pseudotime_range[2], 2), "]\\n")
  
  if (exists("pr_deg_ids")) {
    cat("- Pseudotime DEGs:", length(pr_deg_ids), "\\n")
  }
}

cat("\\nNext steps:\\n")
cat("1. Review trajectory structure and branch points\\n")
cat("2. Validate developmental progression with markers\\n")
cat("3. Perform branch-specific analysis if needed\\n")
cat("4. Integrate with RNA velocity analysis\\n")

cat("\\nUsage examples:\\n")
cat("Rscript monocle3_pseudotime_analysis.R <input_seurat.rds>\\n")