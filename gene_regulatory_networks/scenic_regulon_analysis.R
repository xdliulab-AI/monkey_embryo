#!/usr/bin/env Rscript

# ============================================================================
# SCENIC Regulon Analysis and Activity Scoring
# ============================================================================
# 
# This script performs the complete SCENIC regulon analysis pipeline including
# co-expression module identification, regulon creation, cell scoring, and
# activity binarization.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(SCENIC)
  library(RcisTarget)
  library(doParallel)
  library(dplyr)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(pheatmap)
  library(circlize)
  library(tools)
})

# Parameters
N_CORES <- 1                         # Number of cores for SCENIC steps (memory intensive)
ORGANISM <- "hgnc"                   # "hgnc" for human, "mgi" for mouse
HEATMAP_HEIGHT_FACTOR <- 0.5         # Height scaling factor for heatmaps
HEATMAP_WIDTH <- 12                  # Base width for heatmaps

# Input/Output paths
INT_DIR <- "../output/int"           # Directory with intermediate files
OUTPUT_DIR <- "../output"            # Output directory for results
PLOTS_DIR <- "../plots"              # Directory for plots

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  INT_DIR <- args[1]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to check file existence
check_file_exists <- function(file_path, description) {
  if (!file.exists(file_path)) {
    stop(paste("Required file not found:", description, "-", file_path))
  }
  cat("Found", description, ":", file_path, "\\n")
}

# Function to create comprehensive heatmap
create_regulon_heatmap <- function(activity_matrix, title, filename, 
                                   cell_type_colors = NULL, width = 12, height = 8) {
  
  # Prepare annotation for columns (cell types)
  if (!is.null(cell_type_colors)) {
    col_annotation <- HeatmapAnnotation(
      CellType = colnames(activity_matrix),
      col = list(CellType = cell_type_colors),
      annotation_name_gp = gpar(fontsize = 10)
    )
  } else {
    col_annotation <- NULL
  }
  
  # Create color function for heatmap
  col_fun <- colorRamp2(
    c(min(activity_matrix, na.rm = TRUE), 
      mean(activity_matrix, na.rm = TRUE), 
      max(activity_matrix, na.rm = TRUE)),
    c("blue", "white", "red")
  )
  
  # Create heatmap
  ht <- Heatmap(
    activity_matrix,
    name = "Activity",
    col = col_fun,
    top_annotation = col_annotation,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(
      title = "Scaled\\nActivity",
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  # Save heatmap
  pdf(filename, width = width, height = height)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  return(ht)
}

# Function to analyze regulon targets
analyze_regulon_targets <- function(scenic_options, output_dir) {
  
  # Load regulon target information
  regulon_targets_file <- getIntName(scenic_options, "regulonTargetsInfo")
  
  if (file.exists(regulon_targets_file)) {
    regulon_targets <- loadInt(scenic_options, "regulonTargetsInfo")
    
    # Calculate regulon statistics
    regulon_stats <- regulon_targets %>%
      group_by(TF) %>%
      summarise(
        n_targets = n(),
        mean_importance = mean(importance, na.rm = TRUE),
        median_importance = median(importance, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(desc(n_targets))
    
    # Save regulon statistics
    write.csv(regulon_stats, file.path(output_dir, "regulon_statistics.csv"), 
              row.names = FALSE)
    
    cat("Regulon statistics saved to regulon_statistics.csv\\n")
    return(regulon_stats)
  } else {
    cat("Warning: Regulon targets file not found\\n")
    return(NULL)
  }
}

# ============================================================================
# Load Required Files
# ============================================================================

cat("Starting SCENIC regulon analysis...\\n")
cat("Working directory:", getwd(), "\\n")
cat("Intermediate directory:", INT_DIR, "\\n\\n")

# Check for required files
scenic_options_file <- file.path(INT_DIR, "scenicOptions.Rds")
adjacencies_file <- file.path(INT_DIR, "adj.tsv")
expr_matrix_file <- file.path(INT_DIR, "exprMat_filtered.rds")
cell_info_file <- file.path(INT_DIR, "cellInfo.Rds")
col_vars_file <- file.path(INT_DIR, "colVars.Rds")

# Validate file existence
check_file_exists(scenic_options_file, "SCENIC options")
check_file_exists(adjacencies_file, "Adjacency matrix")
check_file_exists(expr_matrix_file, "Filtered expression matrix")
check_file_exists(cell_info_file, "Cell information")
check_file_exists(col_vars_file, "Color variables")

# ============================================================================
# Load Data and Initialize
# ============================================================================

cat("Loading SCENIC options and data...\\n")

# Load SCENIC options
scenic_options <- readRDS(scenic_options_file)
cat("SCENIC options loaded successfully\\n")

# Load adjacency matrix (GRNBoost output)
cat("Reading GRNBoost/GENIE3 output...\\n")
grn_output <- read.delim(adjacencies_file, header = TRUE)
colnames(grn_output) <- c("TF", "Target", "weight")

# Save as SCENIC intermediate file
genie3_file <- file.path(INT_DIR, "1.4_GENIE3_linkList.Rds")
saveRDS(grn_output, file = genie3_file)
cat("GRNBoost output processed and saved\\n")

# Load expression matrix
cat("Loading filtered expression matrix...\\n")
expr_matrix_filtered <- readRDS(expr_matrix_file)
cat("Expression matrix loaded:", dim(expr_matrix_filtered), "\\n")

# Load cell information and colors
cell_info <- readRDS(cell_info_file)
col_vars <- readRDS(col_vars_file)
cat("Cell information and colors loaded\\n")

# Set number of cores
scenic_options@settings$nCores <- N_CORES
cat("Number of cores set to:", N_CORES, "\\n\\n")

# ============================================================================
# SCENIC Step 1: Gene Correlation
# ============================================================================

cat("========================================\\n")
cat("SCENIC Step 1: Calculating gene correlations\\n")
cat("========================================\\n")

runCorrelation(expr_matrix_filtered, scenic_options)
cat("Gene correlations calculated successfully\\n\\n")

# ============================================================================
# SCENIC Step 2: Co-expression Network to Modules
# ============================================================================

cat("========================================\\n")
cat("SCENIC Step 2: Building co-expression networks and identifying modules\\n")
cat("========================================\\n")

runSCENIC_1_coexNetwork2modules(scenic_options)
cat("Co-expression networks and modules identified\\n\\n")

# ============================================================================
# SCENIC Step 3: Create Regulons
# ============================================================================

cat("========================================\\n")
cat("SCENIC Step 3: Creating regulons\\n")
cat("========================================\\n")

# Load motif annotations based on organism
if (ORGANISM == "hgnc") {
  data("motifAnnotations_hgnc_v9")
  motif_annotations <- motifAnnotations_hgnc_v9
} else {
  data("motifAnnotations_mgi_v9")
  motif_annotations <- motifAnnotations_mgi_v9
}

runSCENIC_2_createRegulons(scenic_options)
cat("Regulons created successfully\\n\\n")

# ============================================================================
# SCENIC Step 4: Score Cells (AUCell)
# ============================================================================

cat("========================================\\n")
cat("SCENIC Step 4: Scoring cells based on regulon activity\\n")
cat("========================================\\n")

# Log-transform expression matrix for AUCell
expr_matrix_log <- log2(expr_matrix_filtered + 1)
cat("Expression matrix log-transformed\\n")

runSCENIC_3_scoreCells(scenic_options, expr_matrix_log)
cat("Cell scoring completed\\n\\n")

# ============================================================================
# SCENIC Step 5: Binarize Regulon Activity
# ============================================================================

cat("========================================\\n")
cat("SCENIC Step 5: Binarizing regulon activity\\n")
cat("========================================\\n")

runSCENIC_4_aucell_binarize(scenic_options)
cat("Regulon activity binarized\\n\\n")

# ============================================================================
# Load and Process Results
# ============================================================================

cat("Loading and processing SCENIC results...\\n")

# Load regulon activity (AUCell scores)
regulon_auc <- loadInt(scenic_options, "aucell_regulonAUC")
regulon_auc <- regulon_auc[onlyNonDuplicatedExtended(rownames(regulon_auc)), ]
cat("Regulon AUC scores loaded:", dim(regulon_auc), "\\n")

# Load binary regulon activity
binary_activity <- loadInt(scenic_options, "aucell_binary_nonDupl")
cat("Binary regulon activity loaded:", dim(binary_activity), "\\n")

# Ensure cell info has proper rownames
if (is.null(rownames(cell_info))) {
  rownames(cell_info) <- colnames(expr_matrix_filtered)
  cat("Assigned rownames to cell_info\\n")
}

# ============================================================================
# Calculate Regulon Activity by Cell Type
# ============================================================================

cat("Calculating regulon activity by cell type...\\n")

# Calculate mean activity by cell type
regulon_activity_by_celltype <- sapply(
  split(rownames(cell_info), cell_info$CellType),
  function(cells) {
    common_cells <- intersect(cells, colnames(regulon_auc))
    if (length(common_cells) > 0) {
      rowMeans(getAUC(regulon_auc)[, common_cells, drop = FALSE])
    } else {
      rep(0, nrow(regulon_auc))
    }
  }
)

# Scale the activity matrix
regulon_activity_scaled <- t(scale(t(regulon_activity_by_celltype), 
                                  center = TRUE, scale = TRUE))

cat("Regulon activity matrix calculated and scaled\\n")
cat("Matrix dimensions:", dim(regulon_activity_scaled), "\\n")

# ============================================================================
# Binary Activity Analysis
# ============================================================================

cat("Analyzing binary regulon activity...\\n")

# Calculate binary activity by cell type
binary_activity_by_celltype <- sapply(
  split(rownames(cell_info), cell_info$CellType),
  function(cells) {
    common_cells <- intersect(cells, colnames(binary_activity))
    if (length(common_cells) > 0) {
      rowMeans(binary_activity[, common_cells, drop = FALSE])
    } else {
      rep(0, nrow(binary_activity))
    }
  }
)

cat("Binary activity matrix calculated\\n")
cat("Matrix dimensions:", dim(binary_activity_by_celltype), "\\n")

# ============================================================================
# Create Visualizations
# ============================================================================

cat("Creating visualizations...\\n")

# Extract cell type colors
cell_type_colors <- col_vars$CellType

# 1. Main regulon activity heatmap (scaled)
main_heatmap_file <- file.path(PLOTS_DIR, "SCENIC_regulonActivity_byCellType_Scaled.pdf")
n_celltypes <- ncol(regulon_activity_scaled)
heatmap_height <- 8 + (n_celltypes * HEATMAP_HEIGHT_FACTOR)

create_regulon_heatmap(
  regulon_activity_scaled,
  title = "Regulon Activity by Cell Type (Scaled)",
  filename = main_heatmap_file,
  cell_type_colors = cell_type_colors,
  width = HEATMAP_WIDTH,
  height = heatmap_height
)

cat("Main regulon activity heatmap saved:", main_heatmap_file, "\\n")

# 2. Binary regulon activity heatmap
binary_heatmap_file <- file.path(PLOTS_DIR, "SCENIC_binaryActivity_byCellType.pdf")

create_regulon_heatmap(
  binary_activity_by_celltype,
  title = "Binary Regulon Activity by Cell Type",
  filename = binary_heatmap_file,
  cell_type_colors = cell_type_colors,
  width = HEATMAP_WIDTH,
  height = heatmap_height
)

cat("Binary regulon activity heatmap saved:", binary_heatmap_file, "\\n")

# 3. Top regulons heatmap (select top 50 most variable)
if (nrow(regulon_activity_scaled) > 50) {
  regulon_variance <- apply(regulon_activity_scaled, 1, var, na.rm = TRUE)
  top_regulons <- names(sort(regulon_variance, decreasing = TRUE)[1:50])
  
  top_heatmap_file <- file.path(PLOTS_DIR, "SCENIC_topRegulons_byCellType.pdf")
  
  create_regulon_heatmap(
    regulon_activity_scaled[top_regulons, ],
    title = "Top 50 Variable Regulons by Cell Type",
    filename = top_heatmap_file,
    cell_type_colors = cell_type_colors,
    width = HEATMAP_WIDTH,
    height = 12
  )
  
  cat("Top regulons heatmap saved:", top_heatmap_file, "\\n")
}

# ============================================================================
# Regulon Analysis and Statistics
# ============================================================================

cat("Analyzing regulon statistics...\\n")

# Analyze regulon targets
regulon_stats <- analyze_regulon_targets(scenic_options, OUTPUT_DIR)

# Calculate cell type specific regulon activity
cell_type_specific_regulons <- data.frame(
  regulon = rownames(regulon_activity_scaled),
  max_activity = apply(regulon_activity_scaled, 1, max, na.rm = TRUE),
  max_celltype = colnames(regulon_activity_scaled)[apply(regulon_activity_scaled, 1, which.max)],
  activity_range = apply(regulon_activity_scaled, 1, function(x) diff(range(x, na.rm = TRUE))),
  stringsAsFactors = FALSE
)

# Sort by activity range (most specific first)
cell_type_specific_regulons <- cell_type_specific_regulons[
  order(cell_type_specific_regulons$activity_range, decreasing = TRUE), 
]

# Save results
write.csv(cell_type_specific_regulons, 
          file.path(OUTPUT_DIR, "cell_type_specific_regulons.csv"), 
          row.names = FALSE)

cat("Cell type specific regulon analysis saved\\n")

# ============================================================================
# Save Results
# ============================================================================

cat("Saving analysis results...\\n")

# Save regulon activity matrices
saveRDS(regulon_activity_by_celltype, 
        file.path(OUTPUT_DIR, "regulon_activity_by_celltype.rds"))
saveRDS(regulon_activity_scaled, 
        file.path(OUTPUT_DIR, "regulon_activity_scaled.rds"))
saveRDS(binary_activity_by_celltype, 
        file.path(OUTPUT_DIR, "binary_activity_by_celltype.rds"))

# Save individual regulon AUC scores
saveRDS(regulon_auc, file.path(OUTPUT_DIR, "regulon_auc_scores.rds"))
saveRDS(binary_activity, file.path(OUTPUT_DIR, "binary_regulon_activity.rds"))

cat("Result matrices saved to output directory\\n")

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating SCENIC analysis report...\\n")

# Calculate summary statistics
n_cells <- ncol(expr_matrix_filtered)
n_genes <- nrow(expr_matrix_filtered)
n_regulons <- nrow(regulon_activity_by_celltype)
n_celltypes <- ncol(regulon_activity_by_celltype)

# Get top active regulons per cell type
top_regulons_per_celltype <- list()
for (ct in colnames(regulon_activity_scaled)) {
  ct_scores <- regulon_activity_scaled[, ct]
  top_5 <- names(sort(ct_scores, decreasing = TRUE)[1:5])
  top_regulons_per_celltype[[ct]] <- top_5
}

# Create report
report_text <- paste0(
  "SCENIC Regulon Analysis Report\\n",
  "=============================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Organism: ", ORGANISM, "\\n\\n",
  "Input Data:\\n",
  "  - Cells: ", n_cells, "\\n",
  "  - Genes: ", n_genes, "\\n",
  "  - Cell types: ", n_celltypes, "\\n\\n",
  "SCENIC Results:\\n",
  "  - Regulons identified: ", n_regulons, "\\n",
  "  - Processing cores: ", N_CORES, "\\n\\n",
  "Analysis Completed:\\n",
  "  1. ✓ Gene correlation calculation\\n",
  "  2. ✓ Co-expression network construction\\n",
  "  3. ✓ Regulon creation with motif analysis\\n",
  "  4. ✓ Cell scoring with AUCell\\n",
  "  5. ✓ Binary activity determination\\n\\n",
  "Output Files:\\n",
  "  - Regulon activity matrix: regulon_activity_by_celltype.rds\\n",
  "  - Scaled activity matrix: regulon_activity_scaled.rds\\n",
  "  - Binary activity matrix: binary_activity_by_celltype.rds\\n",
  "  - Cell type specific regulons: cell_type_specific_regulons.csv\\n"
)

if (!is.null(regulon_stats)) {
  report_text <- paste0(report_text,
    "  - Regulon statistics: regulon_statistics.csv\\n"
  )
}

report_text <- paste0(report_text,
  "\\nVisualization Files:\\n",
  "  - Main heatmap: SCENIC_regulonActivity_byCellType_Scaled.pdf\\n",
  "  - Binary heatmap: SCENIC_binaryActivity_byCellType.pdf\\n"
)

if (nrow(regulon_activity_scaled) > 50) {
  report_text <- paste0(report_text,
    "  - Top regulons: SCENIC_topRegulons_byCellType.pdf\\n"
  )
}

# Add top regulons per cell type
report_text <- paste0(report_text, "\\nTop Active Regulons by Cell Type:\\n")
for (ct in names(top_regulons_per_celltype)) {
  top_regs <- paste(top_regulons_per_celltype[[ct]], collapse = ", ")
  report_text <- paste0(report_text, "  ", ct, ": ", top_regs, "\\n")
}

report_text <- paste0(report_text,
  "\\nNext Steps:\\n",
  "1. Analyze cell type-specific regulons\\n",
  "2. Investigate developmental trajectories\\n",
  "3. Compare regulon activity across conditions\\n",
  "4. Validate key regulons experimentally\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "scenic_regulon_analysis_report.txt"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_scenic_regulon_analysis.txt"))

# ============================================================================
# Final Summary
# ============================================================================

cat("\\n========================================\\n")
cat("SCENIC Regulon Analysis Completed Successfully!\\n")
cat("========================================\\n")
cat("Results summary:\\n")
cat("- Regulons identified:", n_regulons, "\\n")
cat("- Cell types analyzed:", n_celltypes, "\\n")
cat("- Output files saved in:", OUTPUT_DIR, "\\n")
cat("- Plots saved in:", PLOTS_DIR, "\\n")

if (!is.null(regulon_stats)) {
  top_regulon <- regulon_stats$TF[1]
  top_targets <- regulon_stats$n_targets[1]
  cat("- Top regulon:", top_regulon, "(", top_targets, "targets)\\n")
}

cat("\\nMain outputs:\\n")
cat("1. regulon_activity_scaled.rds - Scaled activity matrix\\n")
cat("2. SCENIC_regulonActivity_byCellType_Scaled.pdf - Main heatmap\\n")
cat("3. cell_type_specific_regulons.csv - Regulon specificity analysis\\n")

cat("\\nUsage for next steps:\\n")
cat("Rscript scenic_visualization.R\\n")
cat("Rscript scenic_trajectory_analysis.R\\n")