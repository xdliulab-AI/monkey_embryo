#!/usr/bin/env Rscript

# ============================================================================
# Batch Functional Analysis for Multiple Datasets
# ============================================================================
# 
# This script performs batch functional analysis across multiple datasets,
# clusters, or conditions. It integrates GO enrichment, pathway analysis,
# and custom functional annotation for high-throughput analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(glue)
  library(patchwork)
  library(viridis)
  library(RColorBrewer)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(future)
  library(future.apply)
  library(parallel)
})

# Parameters
JOB_ID <- "batch_functional"
SPECIES <- "human"                    # "human" or "mouse"
BATCH_MODE <- TRUE                    # Process multiple files
P_VALUE_CUTOFF <- 0.01               # P-value cutoff
Q_VALUE_CUTOFF <- 0.01               # Q-value cutoff
MIN_GS_SIZE <- 10                    # Minimum gene set size
MAX_GS_SIZE <- 500                   # Maximum gene set size
N_WORKERS <- 4                       # Parallel workers
TOP_N_TERMS <- 15                    # Top terms per analysis

# Analysis types to perform
ANALYSIS_TYPES <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome")

# Custom pathway categories for filtering
DEVELOPMENTAL_CATEGORIES <- c(
  "embryonic", "development", "morphogenesis", "differentiation",
  "neural", "cardiac", "mesoderm", "endoderm", "ectoderm",
  "signaling", "transcription", "cell cycle", "apoptosis",
  "mesenchymal", "epithelial", "vasculature", "angiogenesis"
)

# Input/Output paths
INPUT_DIR <- "../data"
OUTPUT_DIR <- "../output"
PLOTS_DIR <- "../plots"
DATA_DIR <- "../data"

# File patterns to process
FILE_PATTERNS <- c("*DEG*.rds", "*markers*.rds", "*results*.rds")

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
  SPECIES <- args[2]
}
if (length(args) > 2) {
  JOB_ID <- args[3]
}

# ============================================================================
# Helper Functions
# ============================================================================

# Function to find input files
find_input_files <- function(input_dir, patterns = FILE_PATTERNS) {
  all_files <- c()
  
  for (pattern in patterns) {
    pattern_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    all_files <- c(all_files, pattern_files)
  }
  
  # Remove duplicates and filter for .rds files
  all_files <- unique(all_files)
  all_files <- all_files[grepl("\\.rds$", all_files)]
  
  return(all_files)
}

# Function to extract dataset name from filename
extract_dataset_name <- function(filename) {
  base_name <- basename(filename)
  base_name <- gsub("\\.rds$", "", base_name)
  base_name <- gsub("^[0-9_]*", "", base_name)  # Remove leading numbers
  base_name <- gsub("_(DEG|markers|results).*", "", base_name, ignore.case = TRUE)
  return(base_name)
}

# Function to set up organism info
setup_organism_info <- function(species) {
  if (tolower(species) == "human") {
    return(list(
      org_db = org.Hs.eg.db,
      kegg_organism = "hsa",
      reactome_organism = "human"
    ))
  } else if (tolower(species) == "mouse") {
    return(list(
      org_db = org.Mm.eg.db,
      kegg_organism = "mmu",
      reactome_organism = "mouse"
    ))
  } else {
    stop("Unsupported species. Use 'human' or 'mouse'")
  }
}

# Function to validate and process DEG data
process_deg_data <- function(file_path) {
  
  cat("Processing:", basename(file_path), "\n")
  
  tryCatch({
    data <- readRDS(file_path)
    
    # Validate data structure
    if (!is.data.frame(data)) {
      cat("Warning: Not a data frame, skipping\n")
      return(NULL)
    }
    
    if (!"gene" %in% colnames(data)) {
      cat("Warning: No 'gene' column found, skipping\n")
      return(NULL)
    }
    
    # Remove NA genes
    data <- data[!is.na(data$gene) & data$gene != "", ]
    
    if (nrow(data) == 0) {
      cat("Warning: No valid genes found, skipping\n")
      return(NULL)
    }
    
    cat("Found", nrow(data), "genes\n")
    return(data)
    
  }, error = function(e) {
    cat("Error loading file:", e$message, "\n")
    return(NULL)
  })
}

# Function to convert genes to Entrez IDs
convert_to_entrez <- function(gene_symbols, org_db) {
  
  if (length(gene_symbols) == 0) {
    return(character(0))
  }
  
  tryCatch({
    gene_df <- bitr(gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org_db,
                   drop = TRUE)
    
    return(gene_df$ENTREZID)
    
  }, error = function(e) {
    cat("Gene conversion failed:", e$message, "\n")
    return(character(0))
  })
}

# Function to perform GO enrichment
perform_go_analysis <- function(entrez_ids, org_db, ontology) {
  
  if (length(entrez_ids) < 5) {
    return(NULL)
  }
  
  tryCatch({
    go_result <- enrichGO(
      gene = entrez_ids,
      keyType = "ENTREZID",
      OrgDb = org_db,
      ont = ontology,
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = "fdr",
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      qvalueCutoff = Q_VALUE_CUTOFF,
      readable = TRUE
    )
    
    if (!is.null(go_result) && nrow(go_result) > 0) {
      return(go_result)
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    cat("GO analysis failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to perform KEGG analysis
perform_kegg_analysis <- function(entrez_ids, organism_code) {
  
  if (length(entrez_ids) < 5) {
    return(NULL)
  }
  
  tryCatch({
    kegg_result <- enrichKEGG(
      gene = entrez_ids,
      organism = organism_code,
      keyType = "kegg",
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = "fdr",
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      qvalueCutoff = Q_VALUE_CUTOFF,
      use_internal_data = FALSE
    )
    
    if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
      kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      return(kegg_result)
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    cat("KEGG analysis failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to perform Reactome analysis
perform_reactome_analysis <- function(entrez_ids, organism) {
  
  if (length(entrez_ids) < 5) {
    return(NULL)
  }
  
  tryCatch({
    reactome_result <- enrichPathway(
      gene = entrez_ids,
      organism = organism,
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = "fdr",
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      qvalueCutoff = Q_VALUE_CUTOFF,
      readable = TRUE
    )
    
    if (!is.null(reactome_result) && nrow(reactome_result) > 0) {
      return(reactome_result)
    } else {
      return(NULL)
    }
    
  }, error = function(e) {
    cat("Reactome analysis failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to process single cluster/dataset
process_single_analysis <- function(genes, cluster_name, dataset_name, org_info) {
  
  cat("  Processing cluster:", cluster_name, "\n")
  
  if (length(genes) < 5) {
    cat("  Too few genes, skipping\n")
    return(NULL)
  }
  
  # Convert to Entrez IDs
  entrez_ids <- convert_to_entrez(genes, org_info$org_db)
  
  if (length(entrez_ids) < 5) {
    cat("  Too few valid Entrez IDs, skipping\n")
    return(NULL)
  }
  
  # Perform enrichment analyses
  results <- list()
  
  if ("GO_BP" %in% ANALYSIS_TYPES) {
    go_bp <- perform_go_analysis(entrez_ids, org_info$org_db, "BP")
    if (!is.null(go_bp)) results[["GO_BP"]] <- go_bp
  }
  
  if ("GO_MF" %in% ANALYSIS_TYPES) {
    go_mf <- perform_go_analysis(entrez_ids, org_info$org_db, "MF")
    if (!is.null(go_mf)) results[["GO_MF"]] <- go_mf
  }
  
  if ("GO_CC" %in% ANALYSIS_TYPES) {
    go_cc <- perform_go_analysis(entrez_ids, org_info$org_db, "CC")
    if (!is.null(go_cc)) results[["GO_CC"]] <- go_cc
  }
  
  if ("KEGG" %in% ANALYSIS_TYPES) {
    kegg <- perform_kegg_analysis(entrez_ids, org_info$kegg_organism)
    if (!is.null(kegg)) results[["KEGG"]] <- kegg
  }
  
  if ("Reactome" %in% ANALYSIS_TYPES) {
    reactome <- perform_reactome_analysis(entrez_ids, org_info$reactome_organism)
    if (!is.null(reactome)) results[["Reactome"]] <- reactome
  }
  
  if (length(results) > 0) {
    cat("  Found enriched terms in", length(results), "databases\n")
    return(results)
  } else {
    cat("  No significant enrichment found\n")
    return(NULL)
  }
}

# Function to create summary visualization
create_summary_plot <- function(all_results) {
  
  if (length(all_results) == 0) {
    return(NULL)
  }
  
  # Extract summary statistics
  summary_data <- data.frame(
    Dataset = character(),
    Cluster = character(),
    Database = character(),
    Terms = integer(),
    stringsAsFactors = FALSE
  )
  
  for (dataset_name in names(all_results)) {
    dataset_results <- all_results[[dataset_name]]
    
    for (cluster_name in names(dataset_results)) {
      cluster_results <- dataset_results[[cluster_name]]
      
      for (db_name in names(cluster_results)) {
        n_terms <- nrow(cluster_results[[db_name]])
        
        summary_data <- rbind(summary_data, data.frame(
          Dataset = dataset_name,
          Cluster = cluster_name,
          Database = db_name,
          Terms = n_terms,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(summary_data) == 0) {
    return(NULL)
  }
  
  # Create heatmap plot
  p1 <- ggplot(summary_data, aes(x = paste(Dataset, Cluster, sep = "_"), 
                                y = Database, fill = Terms)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Terms") +
    labs(title = "Functional Enrichment Summary",
         x = "Dataset_Cluster",
         y = "Database") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create bar plot by database
  db_summary <- summary_data %>%
    group_by(Database) %>%
    summarise(Total_Terms = sum(Terms), .groups = "drop")
  
  p2 <- ggplot(db_summary, aes(x = reorder(Database, Total_Terms), y = Total_Terms)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Total Terms by Database",
         x = "Database",
         y = "Total Terms") +
    theme_minimal()
  
  # Combine plots
  combined_plot <- p1 / p2
  
  return(list(plot = combined_plot, data = summary_data))
}

# Function to save batch results
save_batch_results <- function(all_results, output_prefix) {
  
  # Create master Excel workbook
  master_wb <- createWorkbook()
  
  # Summary sheet
  summary_data <- data.frame(
    Dataset = character(),
    Cluster = character(),
    Analysis = character(),
    Terms = integer(),
    Top_Term = character(),
    Top_Pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each dataset
  for (dataset_name in names(all_results)) {
    dataset_results <- all_results[[dataset_name]]
    
    # Create dataset workbook
    dataset_wb <- createWorkbook()
    
    for (cluster_name in names(dataset_results)) {
      cluster_results <- dataset_results[[cluster_name]]
      
      for (analysis_name in names(cluster_results)) {
        result_df <- as.data.frame(cluster_results[[analysis_name]])
        
        # Add to summary
        summary_data <- rbind(summary_data, data.frame(
          Dataset = dataset_name,
          Cluster = cluster_name,
          Analysis = analysis_name,
          Terms = nrow(result_df),
          Top_Term = ifelse(nrow(result_df) > 0, result_df$Description[1], "None"),
          Top_Pvalue = ifelse(nrow(result_df) > 0, result_df$p.adjust[1], 1),
          stringsAsFactors = FALSE
        ))
        
        # Add to dataset workbook
        sheet_name <- paste0(cluster_name, "_", analysis_name)
        sheet_name <- substr(sheet_name, 1, 31)  # Excel limit
        
        addWorksheet(dataset_wb, sheet_name, zoom = 200)
        writeData(dataset_wb, sheet_name, result_df)
        
        # Add to master workbook
        master_sheet_name <- paste0(dataset_name, "_", sheet_name)
        master_sheet_name <- substr(master_sheet_name, 1, 31)
        
        addWorksheet(master_wb, master_sheet_name, zoom = 200)
        writeData(master_wb, master_sheet_name, result_df)
      }
    }
    
    # Save dataset workbook
    dataset_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_", dataset_name, ".xlsx"))
    saveWorkbook(dataset_wb, dataset_file, overwrite = TRUE)
  }
  
  # Add summary to master workbook
  if (nrow(summary_data) > 0) {
    addWorksheet(master_wb, "Summary", zoom = 200)
    writeData(master_wb, "Summary", summary_data)
  }
  
  # Save master workbook
  master_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_master.xlsx"))
  saveWorkbook(master_wb, master_file, overwrite = TRUE)
  
  # Save as RDS
  rds_file <- file.path(DATA_DIR, paste0(output_prefix, "_all_results.rds"))
  saveRDS(all_results, rds_file)
  
  return(list(master = master_file, rds = rds_file, summary = summary_data))
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

cat("Starting batch functional analysis...\n")
cat("Job ID:", JOB_ID, "\n")
cat("Species:", SPECIES, "\n")

# Set up parallel processing
plan(multisession, workers = N_WORKERS)
options(future.globals.maxSize = 3 * 1024^3)

# Set up organism info
org_info <- setup_organism_info(SPECIES)

# Find input files
input_files <- find_input_files(INPUT_DIR, FILE_PATTERNS)
cat("Found", length(input_files), "input files\n")

if (length(input_files) == 0) {
  stop("No input files found in directory: ", INPUT_DIR)
}

# Process each file
all_results <- list()

for (file_path in input_files) {
  dataset_name <- extract_dataset_name(file_path)
  cat("\nProcessing dataset:", dataset_name, "\n")
  
  # Load and validate data
  deg_data <- process_deg_data(file_path)
  
  if (is.null(deg_data)) {
    cat("Skipping", dataset_name, "due to data issues\n")
    next
  }
  
  # Determine analysis strategy
  if ("cluster" %in% colnames(deg_data)) {
    # Cluster-based analysis
    cat("Performing cluster-based analysis\n")
    
    clusters <- unique(deg_data$cluster)
    clusters <- clusters[!is.na(clusters)]
    
    # Limit to manageable number of clusters
    if (length(clusters) > 10) {
      cat("Too many clusters (", length(clusters), "), taking first 10\n")
      clusters <- head(clusters, 10)
    }
    
    dataset_results <- list()
    
    for (cluster in clusters) {
      cluster_genes <- deg_data$gene[deg_data$cluster == cluster]
      cluster_results <- process_single_analysis(cluster_genes, cluster, dataset_name, org_info)
      
      if (!is.null(cluster_results)) {
        dataset_results[[as.character(cluster)]] <- cluster_results
      }
    }
    
    if (length(dataset_results) > 0) {
      all_results[[dataset_name]] <- dataset_results
    }
    
  } else {
    # Single gene set analysis
    cat("Performing single gene set analysis\n")
    
    all_genes <- deg_data$gene
    single_results <- process_single_analysis(all_genes, "all", dataset_name, org_info)
    
    if (!is.null(single_results)) {
      all_results[[dataset_name]] <- list("all" = single_results)
    }
  }
}

cat("\nSuccessfully analyzed", length(all_results), "datasets\n")

# ============================================================================
# Create Visualizations and Save Results
# ============================================================================

if (length(all_results) > 0) {
  
  cat("Creating summary visualizations...\n")
  
  output_prefix <- paste0("batch_functional_", SPECIES, "_", JOB_ID)
  
  # Create summary plot
  summary_result <- create_summary_plot(all_results)
  
  if (!is.null(summary_result)) {
    summary_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_summary.pdf"))
    ggsave(summary_file, summary_result$plot, width = 14, height = 10, dpi = 300)
  }
  
  # Save all results
  cat("Saving results...\n")
  saved_files <- save_batch_results(all_results, output_prefix)
  
  cat("Files saved:\n")
  cat("  Master Excel:", basename(saved_files$master), "\n")
  cat("  RDS file:", basename(saved_files$rds), "\n")
  
} else {
  cat("No results to save\n")
  saved_files <- NULL
  summary_result <- NULL
}

# ============================================================================
# Generate Report
# ============================================================================

cat("Generating comprehensive report...\n")

n_datasets <- length(all_results)
total_analyses <- sum(sapply(all_results, function(x) sum(sapply(x, length))))

report_text <- paste0(
  "Batch Functional Analysis Report\n",
  "===============================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Job ID: ", JOB_ID, "\n",
  "Species: ", SPECIES, "\n",
  "Input Directory: ", INPUT_DIR, "\n",
  "Files Processed: ", length(input_files), "\n",
  "Successful Datasets: ", n_datasets, "\n\n",
  "Parameters:\n",
  "  - P-value cutoff: ", P_VALUE_CUTOFF, "\n",
  "  - Q-value cutoff: ", Q_VALUE_CUTOFF, "\n",
  "  - Min gene set size: ", MIN_GS_SIZE, "\n",
  "  - Max gene set size: ", MAX_GS_SIZE, "\n",
  "  - Analysis types: ", paste(ANALYSIS_TYPES, collapse = ", "), "\n\n",
  "Results Summary:\n",
  "  - Total analyses performed: ", total_analyses, "\n"
)

if (!is.null(summary_result) && !is.null(summary_result$data)) {
  summary_stats <- summary_result$data %>%
    group_by(Database) %>%
    summarise(
      Total_Terms = sum(Terms),
      Avg_Terms = round(mean(Terms), 1),
      .groups = "drop"
    )
  
  report_text <- paste0(report_text, "\nEnrichment by Database:\n")
  
  for (i in 1:nrow(summary_stats)) {
    db_info <- summary_stats[i, ]
    report_text <- paste0(report_text,
      "  - ", db_info$Database, ": ", db_info$Total_Terms, 
      " total terms (avg: ", db_info$Avg_Terms, " per analysis)\n"
    )
  }
}

if (n_datasets > 0) {
  report_text <- paste0(report_text, "\nDatasets Analyzed:\n")
  
  for (dataset_name in names(all_results)) {
    n_clusters <- length(all_results[[dataset_name]])
    report_text <- paste0(report_text,
      "  - ", dataset_name, ": ", n_clusters, " clusters/groups\n"
    )
  }
}

if (!is.null(saved_files)) {
  report_text <- paste0(report_text,
    "\nOutput Files:\n",
    "  - Master Excel: ", basename(saved_files$master), "\n",
    "  - RDS file: ", basename(saved_files$rds), "\n",
    "  - Summary plot: ", paste0(output_prefix, "_summary.pdf"), "\n"
  )
}

report_text <- paste0(report_text,
  "\nRecommendations:\n",
  "1. Review summary heatmap for overall patterns\n",
  "2. Focus on consistently enriched pathways across datasets\n",
  "3. Investigate dataset-specific enrichments\n",
  "4. Validate key functional categories experimentally\n",
  "5. Consider pathway crosstalk and interactions\n\n",
  "Usage Examples:\n",
  "Rscript functional_analysis_batch.R <input_dir> <species> <job_id>\n",
  "Rscript functional_analysis_batch.R data/ human batch_v1\n"
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
cat("Batch functional analysis completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\nFinal Summary:\n")
cat("- Input files processed:", length(input_files), "\n")
cat("- Successful datasets:", n_datasets, "\n")
cat("- Total analyses:", total_analyses, "\n")
cat("- Analysis types:", paste(ANALYSIS_TYPES, collapse = ", "), "\n")

if (!is.null(summary_result) && !is.null(summary_result$data)) {
  most_enriched_db <- summary_result$data %>%
    group_by(Database) %>%
    summarise(Total = sum(Terms), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = 1)
  
  cat("- Most enriched database:", most_enriched_db$Database, 
      "(", most_enriched_db$Total, "terms )\n")
}

cat("\nNext steps:\n")
cat("1. Review batch analysis results and patterns\n")
cat("2. Identify conserved functional themes\n")
cat("3. Investigate dataset-specific enrichments\n")
cat("4. Perform comparative functional analysis\n")

cat("\nAnalysis completed successfully!\n")