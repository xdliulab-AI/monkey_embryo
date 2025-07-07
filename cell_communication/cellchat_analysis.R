#!/usr/bin/env Rscript

# ============================================================================
# CellChat Cell-Cell Communication Analysis
# ============================================================================
# 
# This script performs comprehensive cell-cell communication analysis using
# CellChat to identify ligand-receptor interactions, signaling pathways,
# and communication networks in spatial transcriptomics data.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(gridExtra)
  library(writexl)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(viridis)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
})

# Parameters
SPECIES <- "human"                    # "human" or "mouse"
SIGNALING_TYPE <- "Secreted Signaling"  # "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"
N_WORKERS <- 4                       # Number of parallel workers
TRIM_VALUE <- 0.1                    # Trim value for truncated mean
INTERACTION_RANGE <- 250             # Interaction range for distant signaling
CONTACT_RANGE <- 100                 # Contact range for cell-cell contact
MIN_CELLS <- 10                      # Minimum cells for communication
SPATIAL_RATIO <- 0.5                 # Spatial ratio for spatial factors
SPATIAL_TOLERANCE <- 12.5            # Spatial tolerance

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
if (length(args) > 2) {
  SIGNALING_TYPE <- args[3]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to prepare spatial factors
prepare_spatial_factors <- function(sample_ids, ratio = 0.5, tolerance = 12.5) {
  spatial_factors <- data.frame(
    ratio = rep(ratio, length(sample_ids)),
    tol = rep(tolerance, length(sample_ids))
  )
  rownames(spatial_factors) <- sample_ids
  return(spatial_factors)
}

# Function to create network visualizations
create_network_plots <- function(cellchat_obj, output_prefix, plot_type = "circle") {
  
  # Get group sizes
  group_size <- as.numeric(table(cellchat_obj@idents))
  
  # Create count and weight network plots
  plots_created <- c()
  
  tryCatch({
    # Circle plot for interaction counts
    pdf(file.path(PLOTS_DIR, paste0(output_prefix, "_network_count.pdf")), 
        width = 8, height = 8)
    netVisual_circle(cellchat_obj@net$count, 
                     vertex.weight = rowSums(cellchat_obj@net$count), 
                     weight.scale = TRUE, 
                     label.edge = FALSE, 
                     title.name = "Number of interactions")
    dev.off()
    plots_created <- c(plots_created, "network_count")
    
    # Circle plot for interaction weights
    pdf(file.path(PLOTS_DIR, paste0(output_prefix, "_network_weight.pdf")), 
        width = 8, height = 8)
    netVisual_circle(cellchat_obj@net$weight, 
                     vertex.weight = rowSums(cellchat_obj@net$weight), 
                     weight.scale = TRUE, 
                     label.edge = FALSE, 
                     title.name = "Interaction weights/strength")
    dev.off()
    plots_created <- c(plots_created, "network_weight")
    
  }, error = function(e) {
    cat("Error creating network plots:", e$message, "\\n")
  })
  
  return(plots_created)
}

# Function to create pathway plots
create_pathway_plots <- function(cellchat_obj, output_prefix, max_pathways = 20) {
  
  pathways <- cellchat_obj@netP$pathways
  
  if (length(pathways) == 0) {
    cat("No pathways found for visualization\\n")
    return(NULL)
  }
  
  # Limit number of pathways to avoid huge files
  if (length(pathways) > max_pathways) {
    cat("Too many pathways (", length(pathways), "). Plotting first", max_pathways, "\\n")
    pathways <- pathways[1:max_pathways]
  }
  
  tryCatch({
    # Create PDF with all pathway plots
    pdf(file.path(PLOTS_DIR, paste0(output_prefix, "_pathways.pdf")), 
        width = 10, height = 10)
    
    for (pathway in pathways) {
      tryCatch({
        netVisual_aggregate(cellchat_obj, signaling = pathway, layout = "circle")
        title(main = pathway, cex.main = 1.5)
      }, error = function(e) {
        cat("Error plotting pathway", pathway, ":", e$message, "\\n")
      })
    }
    
    dev.off()
    
    cat("Created pathway plots for", length(pathways), "pathways\\n")
    return("pathways")
    
  }, error = function(e) {
    cat("Error creating pathway plots:", e$message, "\\n")
    return(NULL)
  })
}

# Function to create heatmap visualizations
create_communication_heatmaps <- function(cellchat_obj, output_prefix) {
  
  plots_created <- c()
  
  tryCatch({
    # Communication heatmap
    pdf(file.path(PLOTS_DIR, paste0(output_prefix, "_heatmap_count.pdf")), 
        width = 10, height = 8)
    netVisual_heatmap(cellchat_obj, measure = "count")
    dev.off()
    plots_created <- c(plots_created, "heatmap_count")
    
    # Weight heatmap
    pdf(file.path(PLOTS_DIR, paste0(output_prefix, "_heatmap_weight.pdf")), 
        width = 10, height = 8)
    netVisual_heatmap(cellchat_obj, measure = "weight")
    dev.off()
    plots_created <- c(plots_created, "heatmap_weight")
    
  }, error = function(e) {
    cat("Error creating heatmaps:", e$message, "\\n")
  })
  
  return(plots_created)
}

# Function to create bubble plots
create_bubble_plots <- function(cellchat_obj, output_prefix, top_pathways = 10) {
  
  plots_created <- c()
  
  tryCatch({
    # Get top pathways by communication strength
    pathway_strength <- cellchat_obj@netP$pathways
    
    if (length(pathway_strength) > top_pathways) {
      # Calculate pathway communication probabilities
      pathway_probs <- sapply(pathway_strength, function(p) {
        tryCatch({
          prob_sum <- sum(cellchat_obj@netP$prob[,,p], na.rm = TRUE)
          return(prob_sum)
        }, error = function(e) {
          return(0)
        })
      })
      
      # Select top pathways
      top_pathway_names <- names(sort(pathway_probs, decreasing = TRUE))[1:top_pathways]
    } else {
      top_pathway_names <- pathway_strength
    }
    
    # Create bubble plot
    pdf(file.path(PLOTS_DIR, paste0(output_prefix, "_bubble_plot.pdf")), 
        width = 12, height = 8)
    
    netVisual_bubble(cellchat_obj, sources.use = 1:length(levels(cellchat_obj@idents)), 
                     targets.use = 1:length(levels(cellchat_obj@idents)), 
                     signaling = top_pathway_names, remove.isolate = FALSE)
    
    dev.off()
    plots_created <- c(plots_created, "bubble_plot")
    
  }, error = function(e) {
    cat("Error creating bubble plots:", e$message, "\\n")
  })
  
  return(plots_created)
}

# Function to analyze communication patterns
analyze_communication_patterns <- function(cellchat_obj) {
  
  # Extract communication networks
  net_count <- cellchat_obj@net$count
  net_weight <- cellchat_obj@net$weight
  
  # Calculate basic statistics
  total_interactions <- sum(net_count)
  total_strength <- sum(net_weight)
  
  # Identify top communicating cell types
  outgoing_strength <- rowSums(net_weight)
  incoming_strength <- colSums(net_weight)
  
  communication_stats <- data.frame(
    cell_type = names(outgoing_strength),
    outgoing_strength = outgoing_strength,
    incoming_strength = incoming_strength,
    total_strength = outgoing_strength + incoming_strength,
    outgoing_count = rowSums(net_count),
    incoming_count = colSums(net_count),
    stringsAsFactors = FALSE
  )
  
  communication_stats <- communication_stats[order(communication_stats$total_strength, 
                                                  decreasing = TRUE), ]
  
  # Calculate network metrics
  network_metrics <- list(
    total_interactions = total_interactions,
    total_strength = total_strength,
    avg_interactions_per_pair = total_interactions / (nrow(net_count) * ncol(net_count)),
    top_sender = communication_stats$cell_type[which.max(communication_stats$outgoing_strength)],
    top_receiver = communication_stats$cell_type[which.max(communication_stats$incoming_strength)],
    top_communicator = communication_stats$cell_type[which.max(communication_stats$total_strength)]
  )
  
  return(list(
    communication_stats = communication_stats,
    network_metrics = network_metrics
  ))
}

# ============================================================================
# Load and Prepare Data
# ============================================================================

cat("Loading Seurat object for CellChat analysis...\\n")

if (!file.exists(INPUT_SEURAT)) {
  stop("Input Seurat object not found: ", INPUT_SEURAT)
}

sobj <- readRDS(INPUT_SEURAT)
cat("Loaded Seurat object with", ncol(sobj), "cells and", nrow(sobj), "features\\n")

# Determine cell type column
cell_type_col <- "celltype"
if (!"celltype" %in% colnames(sobj@meta.data)) {
  if ("preliminary_celltype" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "preliminary_celltype"
  } else if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
    cell_type_col <- "seurat_clusters"
  } else {
    stop("No suitable cell type annotation found")
  }
}

cat("Using cell type column:", cell_type_col, "\\n")

# Determine sample column
sample_col <- "sample_id"
if (!"sample_id" %in% colnames(sobj@meta.data)) {
  if ("id" %in% colnames(sobj@meta.data)) {
    sample_col <- "id"
  } else if ("orig.ident" %in% colnames(sobj@meta.data)) {
    sample_col <- "orig.ident"
  } else {
    sample_col <- NULL
    cat("Warning: No sample column found. Using single sample mode.\\n")
  }
}

if (!is.null(sample_col)) {
  cat("Using sample column:", sample_col, "\\n")
}

# ============================================================================
# Prepare CellChat Input Data
# ============================================================================

cat("Preparing CellChat input data...\\n")

# Extract expression data
data_input <- GetAssayData(sobj, slot = "data", assay = "RNA")

# Prepare metadata
if (!is.null(sample_col)) {
  meta <- sobj@meta.data[, c(sample_col, cell_type_col)]
  colnames(meta) <- c("samples", "labels")
} else {
  meta <- data.frame(
    samples = rep("sample1", ncol(sobj)),
    labels = sobj@meta.data[[cell_type_col]]
  )
}

meta$samples <- factor(meta$samples)
meta$labels <- factor(meta$labels)

cat("Cell type distribution:\\n")
print(table(meta$labels))

# Prepare spatial coordinates
spatial_coords <- NULL
coord_source <- "none"

# Try different coordinate sources
if ("Spatial" %in% names(sobj@reductions)) {
  spatial_coords <- Embeddings(sobj, "Spatial")
  coord_source <- "Spatial"
} else if ("Spatial45" %in% names(sobj@reductions)) {
  spatial_coords <- Embeddings(sobj, "Spatial45")
  coord_source <- "Spatial45"
} else if ("SpatialCenter" %in% names(sobj@reductions)) {
  spatial_coords <- Embeddings(sobj, "SpatialCenter")
  coord_source <- "SpatialCenter"
} else if (all(c("x", "y") %in% colnames(sobj@meta.data))) {
  spatial_coords <- as.matrix(sobj@meta.data[, c("x", "y")])
  coord_source <- "metadata"
} else {
  cat("Warning: No spatial coordinates found. Using UMAP coordinates.\\n")
  if ("umap" %in% names(sobj@reductions)) {
    spatial_coords <- Embeddings(sobj, "umap")
    coord_source <- "UMAP"
  }
}

if (!is.null(spatial_coords)) {
  colnames(spatial_coords) <- c("imagerow", "imagecol")
  cat("Using spatial coordinates from:", coord_source, "\\n")
} else {
  stop("No suitable coordinate system found")
}

# Prepare spatial factors
unique_samples <- unique(meta$samples)
spatial_factors <- prepare_spatial_factors(unique_samples, 
                                          ratio = SPATIAL_RATIO, 
                                          tolerance = SPATIAL_TOLERANCE)

cat("Prepared spatial factors for", length(unique_samples), "samples\\n")

# ============================================================================
# Create CellChat Object
# ============================================================================

cat("Creating CellChat object...\\n")

# Create CellChat object
cellchat <- createCellChat(
  object = data_input,
  meta = meta,
  group.by = "labels",
  datatype = "spatial",
  coordinates = spatial_coords,
  spatial.factors = spatial_factors
)

cat("CellChat object created with", ncol(cellchat@data.raw), "cells and", 
    nrow(cellchat@data.raw), "genes\\n")

# ============================================================================
# Set Database and Preprocessing
# ============================================================================

cat("Setting up ligand-receptor database...\\n")

# Load appropriate database
if (SPECIES == "human") {
  CellChatDB <- CellChatDB.human
} else if (SPECIES == "mouse") {
  CellChatDB <- CellChatDB.mouse
} else {
  stop("Unsupported species: ", SPECIES)
}

# Subset database based on signaling type
if (SIGNALING_TYPE %in% c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")) {
  CellChatDB.use <- subsetDB(CellChatDB, search = SIGNALING_TYPE, key = "annotation")
  cat("Using", SIGNALING_TYPE, "database with", nrow(CellChatDB.use$interaction), "interactions\\n")
} else {
  CellChatDB.use <- CellChatDB
  cat("Using full CellChat database with", nrow(CellChatDB.use$interaction), "interactions\\n")
}

# Set database in object
cellchat@DB <- CellChatDB.use

# ============================================================================
# Preprocessing and Computation
# ============================================================================

cat("Preprocessing expression data...\\n")

# Subset data
cellchat <- subsetData(cellchat)

# Set up parallel processing
future::plan("multisession", workers = N_WORKERS)
options(future.globals.maxSize = 3 * 1024^3)  # 3 GB

# Identify over-expressed genes and interactions
cat("Identifying over-expressed genes and interactions...\\n")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities
cat("Computing communication probabilities...\\n")
cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = TRIM_VALUE,
  distance.use = TRUE,
  interaction.range = INTERACTION_RANGE,
  scale.distance = NULL,
  contact.dependent = TRUE,
  contact.range = CONTACT_RANGE
)

# Filter communications
cat("Filtering communications...\\n")
cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)

# Compute pathway communication probabilities
cat("Computing pathway communication probabilities...\\n")
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate networks
cat("Aggregating communication networks...\\n")
cellchat <- aggregateNet(cellchat)

cat("CellChat analysis completed successfully\\n")

# ============================================================================
# Analysis and Visualization
# ============================================================================

cat("Creating visualizations and analyzing results...\\n")

# Create output prefix
output_prefix <- paste0("cellchat_", SPECIES, "_", gsub(" ", "_", SIGNALING_TYPE))

# 1. Network visualizations
network_plots <- create_network_plots(cellchat, output_prefix)

# 2. Pathway visualizations
pathway_plots <- create_pathway_plots(cellchat, output_prefix, max_pathways = 15)

# 3. Heatmap visualizations
heatmap_plots <- create_communication_heatmaps(cellchat, output_prefix)

# 4. Bubble plots
bubble_plots <- create_bubble_plots(cellchat, output_prefix, top_pathways = 10)

# 5. Analyze communication patterns
communication_analysis <- analyze_communication_patterns(cellchat)

# ============================================================================
# Extract and Save Results
# ============================================================================

cat("Extracting and saving results...\\n")

# Extract ligand-receptor interactions
lr_interactions <- subsetCommunication(cellchat)

if (nrow(lr_interactions) > 0) {
  # Save as Excel file
  write_xlsx(lr_interactions, file.path(OUTPUT_DIR, paste0(output_prefix, "_LR_interactions.xlsx")))
  
  # Save as CSV
  write.csv(lr_interactions, file.path(OUTPUT_DIR, paste0(output_prefix, "_LR_interactions.csv")), 
            row.names = FALSE)
  
  cat("Saved", nrow(lr_interactions), "ligand-receptor interactions\\n")
} else {
  cat("Warning: No ligand-receptor interactions found\\n")
}

# Save communication statistics
write.csv(communication_analysis$communication_stats, 
          file.path(OUTPUT_DIR, paste0(output_prefix, "_communication_stats.csv")), 
          row.names = FALSE)

# Save pathway information
if (length(cellchat@netP$pathways) > 0) {
  pathway_info <- data.frame(
    pathway = cellchat@netP$pathways,
    stringsAsFactors = FALSE
  )
  
  # Calculate pathway statistics
  pathway_info$total_interactions <- sapply(pathway_info$pathway, function(p) {
    tryCatch({
      sum(cellchat@netP$count[,,p], na.rm = TRUE)
    }, error = function(e) { 0 })
  })
  
  pathway_info$total_strength <- sapply(pathway_info$pathway, function(p) {
    tryCatch({
      sum(cellchat@netP$prob[,,p], na.rm = TRUE)
    }, error = function(e) { 0 })
  })
  
  pathway_info <- pathway_info[order(pathway_info$total_strength, decreasing = TRUE), ]
  
  write.csv(pathway_info, file.path(OUTPUT_DIR, paste0(output_prefix, "_pathway_info.csv")), 
            row.names = FALSE)
}

# Save CellChat object
saveRDS(cellchat, file.path(OUTPUT_DIR, paste0(output_prefix, "_cellchat_object.rds")))

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("Generating comprehensive report...\\n")

# Calculate summary statistics
n_cells <- ncol(cellchat@data.raw)
n_cell_types <- length(unique(cellchat@idents))
n_interactions <- nrow(lr_interactions)
n_pathways <- length(cellchat@netP$pathways)

report_text <- paste0(
  "CellChat Cell-Cell Communication Analysis Report\\n",
  "==============================================\\n\\n",
  "Analysis Date: ", Sys.Date(), "\\n",
  "Species: ", SPECIES, "\\n",
  "Signaling Type: ", SIGNALING_TYPE, "\\n\\n",
  "Data Summary:\\n",
  "  - Total cells: ", n_cells, "\\n",
  "  - Cell types: ", n_cell_types, "\\n",
  "  - Coordinate source: ", coord_source, "\\n",
  "  - Samples: ", length(unique_samples), "\\n\\n",
  "Analysis Parameters:\\n",
  "  - Interaction range: ", INTERACTION_RANGE, "\\n",
  "  - Contact range: ", CONTACT_RANGE, "\\n",
  "  - Minimum cells: ", MIN_CELLS, "\\n",
  "  - Trim value: ", TRIM_VALUE, "\\n",
  "  - Spatial ratio: ", SPATIAL_RATIO, "\\n",
  "  - Spatial tolerance: ", SPATIAL_TOLERANCE, "\\n\\n",
  "Results Summary:\\n",
  "  - Ligand-receptor interactions: ", n_interactions, "\\n",
  "  - Signaling pathways: ", n_pathways, "\\n"
)

# Add network metrics
if (!is.null(communication_analysis$network_metrics)) {
  nm <- communication_analysis$network_metrics
  report_text <- paste0(report_text,
    "  - Total interaction strength: ", round(nm$total_strength, 2), "\\n",
    "  - Top sender: ", nm$top_sender, "\\n",
    "  - Top receiver: ", nm$top_receiver, "\\n",
    "  - Top communicator: ", nm$top_communicator, "\\n"
  )
}

# Add top pathways
if (exists("pathway_info") && nrow(pathway_info) > 0) {
  report_text <- paste0(report_text,
    "\\nTop 5 Signaling Pathways:\\n"
  )
  
  top_pathways <- head(pathway_info, 5)
  for (i in 1:nrow(top_pathways)) {
    report_text <- paste0(report_text,
      "  ", i, ". ", top_pathways$pathway[i], " (strength: ", 
      round(top_pathways$total_strength[i], 3), ")\\n"
    )
  }
}

report_text <- paste0(report_text,
  "\\nOutput Files:\\n",
  "  - CellChat object: ", output_prefix, "_cellchat_object.rds\\n",
  "  - L-R interactions: ", output_prefix, "_LR_interactions.xlsx\\n",
  "  - Communication stats: ", output_prefix, "_communication_stats.csv\\n"
)

if (exists("pathway_info")) {
  report_text <- paste0(report_text,
    "  - Pathway info: ", output_prefix, "_pathway_info.csv\\n"
  )
}

# Add visualization files
all_plots <- c(network_plots, pathway_plots, heatmap_plots, bubble_plots)
if (length(all_plots) > 0) {
  report_text <- paste0(report_text,
    "\\nVisualization Files:\\n"
  )
  for (plot_type in all_plots) {
    report_text <- paste0(report_text,
      "  - ", plot_type, ": ", output_prefix, "_", plot_type, ".pdf\\n"
    )
  }
}

report_text <- paste0(report_text,
  "\\nRecommendations:\\n",
  "1. Review top communicating cell types for biological relevance\\n",
  "2. Validate key ligand-receptor pairs with expression data\\n",
  "3. Investigate pathway activities in developmental context\\n",
  "4. Compare communication patterns across conditions\\n",
  "5. Perform functional analysis of communication networks\\n"
)

writeLines(report_text, file.path(OUTPUT_DIR, paste0(output_prefix, "_analysis_report.txt")))

# Save session information
writeLines(capture.output(sessionInfo()), 
           file.path(OUTPUT_DIR, paste0(output_prefix, "_sessionInfo.txt")))

# ============================================================================
# Final Summary
# ============================================================================

cat("\\nCellChat analysis completed successfully!\\n")
cat("Results saved in:", OUTPUT_DIR, "\\n")
cat("Plots saved in:", PLOTS_DIR, "\\n")
cat("\\nSummary:\\n")
cat("- Cells analyzed:", n_cells, "\\n")
cat("- Cell types:", n_cell_types, "\\n")
cat("- L-R interactions:", n_interactions, "\\n")
cat("- Signaling pathways:", n_pathways, "\\n")

if (!is.null(communication_analysis$network_metrics)) {
  cat("- Top communicator:", communication_analysis$network_metrics$top_communicator, "\\n")
}

cat("\\nNext steps:\\n")
cat("1. Review communication networks and validate key interactions\\n")
cat("2. Analyze pathway activities in biological context\\n")
cat("3. Compare communication patterns across samples/conditions\\n")
cat("4. Investigate spatial organization of communicating cells\\n")

cat("\\nUsage examples:\\n")
cat("Rscript cellchat_analysis.R <input_seurat.rds> <species> <signaling_type>\\n")
cat("Rscript cellchat_analysis.R data.rds human 'Secreted Signaling'\\n")
cat("Rscript cellchat_analysis.R data.rds mouse 'ECM-Receptor'\\n")