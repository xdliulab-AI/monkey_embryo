#!/usr/bin/env Rscript

# ============================================================================
# SCENIC Data Preparation for Gene Regulatory Network Analysis
# ============================================================================
# 
# This script prepares single-cell RNA-seq data for SCENIC analysis by
# performing gene filtering, cell sampling, and creating the necessary
# input files for gene regulatory network inference.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(dplyr)
  library(SCENIC)
  library(RcisTarget)
  library(SCopeLoomR)
  library(RColorBrewer)
  library(tools)
})

# Parameters
ORGANISM <- "hgnc"                    # "hgnc" for human, "mgi" for mouse
N_CELLS_SAMPLE <- 10000              # Number of cells to sample per dataset
N_CORES <- 12                        # Number of cores for parallel processing
MIN_COUNTS_PER_GENE_PERCENT <- 0.01  # Minimum percentage for gene filtering
MIN_SAMPLES_PERCENT <- 0.01          # Minimum percentage of samples for gene filtering
DATASET_TITLE <- "Embryo_SCENIC_Analysis"

# Input/Output paths
INPUT_SEURAT <- "../data/annotated_seurat_object.rds"
CISTARGGET_DB_DIR <- "../data/cisTarget_databases/hg38"  # Path to cisTarget databases
OUTPUT_DIR <- "../output"
INT_DIR <- file.path(OUTPUT_DIR, "int")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(INT_DIR, showWarnings = FALSE, recursive = TRUE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  INPUT_SEURAT <- args[1]
}
if (length(args) > 1) {
  N_CELLS_SAMPLE <- as.numeric(args[2])
}
if (length(args) > 2) {
  CISTARGGET_DB_DIR <- args[3]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to create color palette for cell types
create_cell_type_colors <- function(cell_types) {
  num_colors <- length(cell_types)
  
  if (num_colors <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, num_colors), "Set1")[1:num_colors]
  } else if (num_colors <= 12) {
    colors <- RColorBrewer::brewer.pal(num_colors, "Set3")
  } else {
    # Use a combination of color palettes for more than 12 colors
    colors <- c(
      RColorBrewer::brewer.pal(8, "Set1"),
      RColorBrewer::brewer.pal(min(num_colors - 8, 12), "Set3")
    )
    if (num_colors > 20) {
      colors <- colorRampPalette(colors)(num_colors)
    }
  }
  
  return(setNames(colors, cell_types))
}

# Function to validate cisTarget databases
validate_cistarget_db <- function(db_dir, organism) {
  if (!dir.exists(db_dir)) {
    stop("cisTarget database directory not found: ", db_dir)
  }
  
  # Check for required database files
  required_files <- c(
    paste0(organism, "__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"),
    paste0(organism, "__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
  )
  
  missing_files <- c()
  for (file in required_files) {
    full_path <- file.path(db_dir, file)
    if (!file.exists(full_path)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    cat("Warning: Missing cisTarget database files:\\n")
    cat(paste(missing_files, collapse = "\\n"), "\\n")
    cat("SCENIC analysis may not work properly without these files.\\n")
  } else {
    cat("All required cisTarget database files found.\\n")
  }
}

# ============================================================================
# Load and Validate Data
# ============================================================================

cat("Starting SCENIC data preparation...\\n")
cat("Dataset title:", DATASET_TITLE, "\\n")
cat("Organism:", ORGANISM, "\\n")
cat("Cells to sample:", N_CELLS_SAMPLE, "\\n")
cat("Number of cores:", N_CORES, "\\n\\n")

# Load Seurat object
cat("Loading Seurat object from:", INPUT_SEURAT, "\\n")
if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded Seurat object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Validate cisTarget databases
validate_cistarget_db(CISTARGGET_DB_DIR, ifelse(ORGANISM == "hgnc", "hg38", "mm10"))

# ============================================================================
# Cell Type and Metadata Preparation
# ============================================================================

cat("\\nPreparing cell type annotations...\\n")

# Determine cell type column
cell_type_col <- "celltype"
if (!"celltype" %in% colnames(sobj@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "preliminary_celltype"
  } else if ("cell_type" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "cell_type"
  } else if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "seurat_clusters"
    cat("Warning: Using cluster IDs as cell types\\n")
  } else {
    stop("No suitable cell type annotation found")
  }
}

cat("Using cell type column:", cell_type_col, "\\n")

# Create minimal metadata for SCENIC
sobj@meta.data <- sobj@meta.data[, c("orig.ident", cell_type_col)]
colnames(sobj@meta.data)[colnames(sobj@meta.data) == cell_type_col] <- "CellType"

# Set identity
Idents(sobj) <- "CellType"

# Remove cells with missing cell type annotations
valid_cells <- !is.na(sobj$CellType) & sobj$CellType != "" & sobj$CellType != "Unknown"
if (sum(!valid_cells) > 0) {
  cat("Removing", sum(!valid_cells), "cells with missing or unknown cell type annotations\\n")
  sobj <- subset(sobj, cells = colnames(sobj)[valid_cells])
}

# Get unique cell types
unique_cell_types <- unique(sobj$CellType)
unique_cell_types <- unique_cell_types[!is.na(unique_cell_types)]
cat("Found", length(unique_cell_types), "unique cell types:\\n")
cat(paste(unique_cell_types, collapse = ", "), "\\n")

# ============================================================================
# Cell Sampling
# ============================================================================

cat("\\nPerforming cell sampling...\\n")

# Check if sampling is needed
total_cells <- ncol(sobj)
if (total_cells > N_CELLS_SAMPLE) {
  cat("Sampling", N_CELLS_SAMPLE, "cells from", total_cells, "total cells\\n")
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Stratified sampling to maintain cell type proportions
  cell_type_counts <- table(sobj$CellType)
  cell_type_props <- cell_type_counts / sum(cell_type_counts)
  
  selected_cells <- c()
  for (ct in names(cell_type_props)) {
    ct_cells <- colnames(sobj)[sobj$CellType == ct]
    n_sample_ct <- round(N_CELLS_SAMPLE * cell_type_props[ct])
    n_sample_ct <- min(n_sample_ct, length(ct_cells))  # Don't exceed available cells
    
    if (n_sample_ct > 0) {
      sampled_ct_cells <- sample(ct_cells, n_sample_ct)
      selected_cells <- c(selected_cells, sampled_ct_cells)
    }
  }
  
  # If we still need more cells, randomly sample from remaining
  if (length(selected_cells) < N_CELLS_SAMPLE) {
    remaining_cells <- setdiff(colnames(sobj), selected_cells)
    additional_needed <- N_CELLS_SAMPLE - length(selected_cells)
    additional_cells <- sample(remaining_cells, min(additional_needed, length(remaining_cells)))
    selected_cells <- c(selected_cells, additional_cells)
  }
  
  # Subset Seurat object
  sobj <- subset(sobj, cells = selected_cells)
  
} else {
  cat("Using all", total_cells, "cells (less than sampling threshold)\\n")
}

cat("Final dataset contains", ncol(sobj), "cells\\n")

# Print final cell type distribution
cat("\\nFinal cell type distribution:\\n")
print(table(sobj$CellType))

# ============================================================================
# Expression Matrix Preparation
# ============================================================================

cat("\\nPreparing expression matrix...\\n")

# Extract count matrix
expr_matrix <- LayerData(sobj, assay = "RNA", layer = "counts")
expr_matrix <- as.matrix(expr_matrix)

cat("Expression matrix dimensions:", dim(expr_matrix), "\\n")

# ============================================================================
# SCENIC Options Initialization
# ============================================================================

cat("\\nInitializing SCENIC options...\\n")

# Load default database names
data(defaultDbNames)

# Set database names for the organism
if (ORGANISM == "hgnc") {
  defaultDbNames$hgnc["500bp"] <- "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
  defaultDbNames$hgnc["10kb"] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
  db_names <- defaultDbNames[["hgnc"]]
} else if (ORGANISM == "mgi") {
  db_names <- defaultDbNames[["mgi"]]
} else {
  stop("Unsupported organism: ", ORGANISM)
}

# Load motif annotations
if (ORGANISM == "hgnc") {
  data("motifAnnotations_hgnc_v9")
  motif_annotations <- motifAnnotations_hgnc_v9
} else {
  data("motifAnnotations_mgi_v9")
  motif_annotations <- motifAnnotations_mgi_v9
}

# Initialize SCENIC options
scenic_options <- initializeScenic(
  org = ORGANISM,
  dbDir = CISTARGGET_DB_DIR,
  dbs = db_names,
  datasetTitle = DATASET_TITLE,
  nCores = N_CORES
)

# Set additional options
scenic_options@settings$verbose <- TRUE
scenic_options@settings$seed <- 123

cat("SCENIC options initialized successfully\\n")

# ============================================================================
# Cell Information and Color Setup
# ============================================================================

cat("\\nPreparing cell information and colors...\\n")

# Create cell information dataframe
cell_info <- sobj@meta.data %>%
  select(orig.ident, CellType) %>%
  mutate(CellType = factor(CellType))

# Drop unused factor levels
cell_info$CellType <- droplevels(cell_info$CellType)

# Create color variables
unique_cell_types_final <- levels(cell_info$CellType)
cell_type_colors <- create_cell_type_colors(unique_cell_types_final)

col_vars <- list(CellType = cell_type_colors)

cat("Created color palette for", length(unique_cell_types_final), "cell types\\n")

# ============================================================================
# Gene Filtering
# ============================================================================

cat("\\nPerforming gene filtering...\\n")

# Calculate filtering thresholds
min_counts_per_gene <- 3 * MIN_COUNTS_PER_GENE_PERCENT * ncol(expr_matrix)
min_samples <- ncol(expr_matrix) * MIN_SAMPLES_PERCENT

cat("Gene filtering parameters:\\n")
cat("  - Minimum counts per gene:", round(min_counts_per_gene, 2), "\\n")
cat("  - Minimum samples per gene:", round(min_samples, 2), "\\n")

# Apply gene filtering
genes_kept <- geneFiltering(
  expr_matrix,
  scenicOptions = scenic_options,
  minCountsPerGene = min_counts_per_gene,
  minSamples = min_samples
)

expr_matrix_filtered <- expr_matrix[genes_kept, ]

cat("Genes kept after filtering:", length(genes_kept), "/", nrow(expr_matrix), 
    "(", round(100 * length(genes_kept) / nrow(expr_matrix), 1), "%)\\n")
cat("Filtered expression matrix dimensions:", dim(expr_matrix_filtered), "\\n")

# ============================================================================
# Save Intermediate Files
# ============================================================================

cat("\\nSaving intermediate files...\\n")

# Save cell information
cell_info_file <- file.path(INT_DIR, "cellInfo.Rds")
saveRDS(cell_info, file = cell_info_file)
cat("Cell information saved to:", cell_info_file, "\\n")

# Save color variables
col_vars_file <- file.path(INT_DIR, "colVars.Rds")
saveRDS(col_vars, file = col_vars_file)
cat("Color variables saved to:", col_vars_file, "\\n")

# Update SCENIC options with file paths
scenic_options@inputDatasetInfo$cellInfo <- cell_info_file
scenic_options@inputDatasetInfo$colVars <- col_vars_file

# Save SCENIC options
scenic_options_file <- file.path(INT_DIR, "scenicOptions.Rds")
saveRDS(scenic_options, file = scenic_options_file)
cat("SCENIC options saved to:", scenic_options_file, "\\n")

# Save filtered expression matrix
expr_matrix_file <- file.path(INT_DIR, "exprMat_filtered.rds")
saveRDS(expr_matrix_filtered, file = expr_matrix_file)
cat("Filtered expression matrix saved to:", expr_matrix_file, "\\n")

# ============================================================================
# Create Loom File
# ============================================================================

cat("\\nCreating loom file for network inference...\\n")

# Build loom file
loom_file <- file.path(INT_DIR, paste0(DATASET_TITLE, "_exprMat_filtered.loom"))

tryCatch({
  loom <- build_loom(
    loom_file,
    dgem = expr_matrix_filtered
  )
  close_loom(loom)
  cat("Loom file created successfully:", loom_file, "\\n")
}, error = function(e) {
  cat("Error creating loom file:", e$message, "\\n")
  cat("This may affect network inference step\\n")
})

# ============================================================================
# Generate Summary Report
# ============================================================================

cat("\\nGenerating data preparation report...\\n")

# Calculate summary statistics
n_final_cells <- ncol(expr_matrix_filtered)
n_final_genes <- nrow(expr_matrix_filtered)
n_cell_types <- length(unique_cell_types_final)

# Create summary report
report_text <- paste0(
  "SCENIC Data Preparation Report\\n",
  "=============================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Dataset Title: ", DATASET_TITLE, "\\n",
  "Organism: ", ORGANISM, "\\n\\n",
  "Input Data:\\n",
  "  - Original cells: ", total_cells, "\\n",
  "  - Original genes: ", nrow(expr_matrix), "\\n",
  "  - Sampling target: ", N_CELLS_SAMPLE, "\\n\\n",
  "Final Data for SCENIC:\\n",
  "  - Cells: ", n_final_cells, "\\n",
  "  - Genes: ", n_final_genes, "\\n",
  "  - Cell types: ", n_cell_types, "\\n\\n",
  "Gene Filtering:\\n",
  "  - Min counts per gene: ", round(min_counts_per_gene, 2), "\\n",
  "  - Min samples: ", round(min_samples, 2), "\\n",
  "  - Genes retained: ", round(100 * n_final_genes / nrow(expr_matrix), 1), "%\\n\\n",
  "Cell Type Distribution:\\n"
)

# Add cell type counts to report
cell_type_table <- table(cell_info$CellType)
for (ct in names(cell_type_table)) {
  report_text <- paste0(report_text,
    "  - ", ct, ": ", cell_type_table[ct], " cells\\n"
  )
}

report_text <- paste0(report_text,
  "\\nOutput Files:\\n",
  "  - Cell info: cellInfo.Rds\\n",
  "  - Color variables: colVars.Rds\\n",
  "  - SCENIC options: scenicOptions.Rds\\n",
  "  - Filtered expression matrix: exprMat_filtered.rds\\n",
  "  - Loom file: ", basename(loom_file), "\\n\\n",
  "Next Steps:\\n",
  "1. Run network inference using GRNBoost2 or GENIE3\\n",
  "2. Execute SCENIC regulon analysis pipeline\\n",
  "3. Analyze regulon activity and cell state transitions\\n\\n",
  "Required for Next Step:\\n",
  "  - cisTarget databases in: ", CISTARGGET_DB_DIR, "\\n",
  "  - Transcription factor list for ", ORGANISM, "\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, "scenic_data_preparation_report.txt"))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, "sessionInfo_scenic_preparation.txt"))

# ============================================================================
# Final Summary
# ============================================================================

cat("\\nSCENIC data preparation completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Intermediate files in:", INT_DIR, "\\n")
cat("\\nSummary:\\n")
cat("- Final dataset:", n_final_cells, "cells x", n_final_genes, "genes\\n")
cat("- Cell types:", n_cell_types, "\\n")
cat("- Gene filtering efficiency:", round(100 * n_final_genes / nrow(expr_matrix), 1), "%\\n")

if (file.exists(loom_file)) {
  cat("- Loom file ready for network inference\\n")
}

cat("\\nNext steps:\\n")
cat("1. Run network inference: scenic_network_inference.py\\n")
cat("2. Generate regulons: scenic_regulon_analysis.R\\n")
cat("3. Analyze results: scenic_visualization.R\\n")

cat("\\nUsage examples:\\n")
cat("Rscript scenic_data_preparation.R <seurat_file> <n_cells> <db_dir>\\n")
cat("Rscript scenic_data_preparation.R data.rds 5000 /path/to/cistargget/\\n")