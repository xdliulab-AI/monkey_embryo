#!/usr/bin/env Rscript

# ============================================================================
# Comprehensive Gene Ontology (GO) Enrichment Analysis
# ============================================================================
# 
# This script performs comprehensive GO enrichment analysis on differentially
# expressed genes from single-cell or spatial transcriptomics data. It includes
# multiple visualization options, pathway filtering, and statistical analysis.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
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
})

# Parameters
JOB_ID <- "go_analysis"
SPECIES <- "human"                    # "human" or "mouse"
ONTOLOGY <- "ALL"                     # "ALL", "BP", "MF", "CC"
P_VALUE_CUTOFF <- 0.01               # P-value cutoff for significance
Q_VALUE_CUTOFF <- 0.01               # Q-value cutoff for significance
MIN_GS_SIZE <- 10                    # Minimum gene set size
MAX_GS_SIZE <- 500                   # Maximum gene set size
P_ADJUST_METHOD <- "fdr"             # P-value adjustment method
LOG2FC_THRESHOLD <- 0.5              # Log2FC threshold for gene filtering
TOP_N_TERMS <- 20                    # Number of top terms to show in plots

# Input/Output paths
INPUT_DEG_FILE <- "../data/deg_results.rds"
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
  INPUT_DEG_FILE <- args[1]
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

# Function to set up organism database
setup_organism_db <- function(species) {
  if (tolower(species) == "human") {
    return(org.Hs.eg.db)
  } else if (tolower(species) == "mouse") {
    return(org.Mm.eg.db)
  } else {
    stop("Unsupported species. Use 'human' or 'mouse'")
  }
}

# Function to convert gene symbols to Entrez IDs
convert_gene_ids <- function(gene_symbols, org_db) {
  
  cat("Converting", length(gene_symbols), "gene symbols to Entrez IDs...\n")
  
  # Remove duplicates and filter out empty/NA genes
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  
  if (length(gene_symbols) == 0) {
    stop("No valid gene symbols provided")
  }
  
  # Convert gene symbols to Entrez IDs
  tryCatch({
    gene_df <- bitr(gene_symbols, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org_db,
                   drop = TRUE)
    
    cat("Successfully converted", nrow(gene_df), "genes\n")
    
    # Check conversion rate
    conversion_rate <- nrow(gene_df) / length(gene_symbols)
    cat("Conversion rate:", round(conversion_rate * 100, 1), "%\n")
    
    if (conversion_rate < 0.5) {
      warning("Low gene ID conversion rate. Check gene symbols and organism database.")
    }
    
    return(gene_df)
    
  }, error = function(e) {
    stop("Gene ID conversion failed: ", e$message)
  })
}

# Function to perform GO enrichment analysis
perform_go_enrichment <- function(entrez_ids, org_db, ontology = ONTOLOGY) {
  
  cat("Performing GO enrichment analysis...\n")
  cat("Genes:", length(entrez_ids), "\n")
  cat("Ontology:", ontology, "\n")
  
  if (length(entrez_ids) == 0) {
    stop("No Entrez IDs provided for enrichment analysis")
  }
  
  tryCatch({
    go_results <- enrichGO(
      gene = entrez_ids,
      keyType = "ENTREZID",
      OrgDb = org_db,
      ont = ontology,
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = P_ADJUST_METHOD,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      qvalueCutoff = Q_VALUE_CUTOFF,
      readable = TRUE
    )
    
    if (!is.null(go_results) && nrow(go_results) > 0) {
      cat("Found", nrow(go_results), "significant GO terms\n")
      return(go_results)
    } else {
      cat("No significant GO terms found\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("GO enrichment analysis failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to filter GO results by custom criteria
filter_go_results <- function(go_results, custom_terms = NULL, min_count = 3) {
  
  if (is.null(go_results) || nrow(go_results) == 0) {
    return(NULL)
  }
  
  go_df <- as.data.frame(go_results)
  
  # Filter by minimum gene count
  go_df <- go_df[go_df$Count >= min_count, ]
  
  # Filter by custom terms if provided
  if (!is.null(custom_terms) && length(custom_terms) > 0) {
    # Create pattern to match any of the custom terms
    pattern <- paste(custom_terms, collapse = "|")
    go_df <- go_df[grepl(pattern, go_df$Description, ignore.case = TRUE), ]
    cat("Filtered to", nrow(go_df), "terms matching custom criteria\n")
  }
  
  return(go_df)
}

# Function to create bar plot
create_go_barplot <- function(go_results, title = "GO Enrichment", top_n = TOP_N_TERMS) {
  
  if (is.null(go_results) || nrow(go_results) == 0) {
    return(NULL)
  }
  
  go_df <- as.data.frame(go_results)
  
  # Select top terms
  if (nrow(go_df) > top_n) {
    go_df <- head(go_df, top_n)
  }
  
  # Create plot
  p <- go_df %>%
    mutate(Description = str_wrap(Description, width = 50)) %>%
    mutate(Description = fct_reorder(Description, Count)) %>%
    ggplot(aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), alpha = 0.8) +
    scale_fill_viridis_c(name = "Adjusted\nP-value", trans = "log10",
                         guide = guide_colorbar(reverse = TRUE)) +
    coord_flip() +
    labs(title = title,
         x = "GO Terms", 
         y = "Gene Count") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 10))
  
  return(p)
}

# Function to create dot plot
create_go_dotplot <- function(go_results, title = "GO Enrichment", top_n = TOP_N_TERMS) {
  
  if (is.null(go_results) || nrow(go_results) == 0) {
    return(NULL)
  }
  
  go_df <- as.data.frame(go_results)
  
  # Select top terms
  if (nrow(go_df) > top_n) {
    go_df <- head(go_df, top_n)
  }
  
  # Calculate gene ratio
  go_df$GeneRatio_numeric <- sapply(go_df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  
  # Create plot
  p <- go_df %>%
    mutate(Description = str_wrap(Description, width = 50)) %>%
    mutate(Description = fct_reorder(Description, Count)) %>%
    ggplot(aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.7) +
    scale_color_viridis_c(name = "Adjusted\nP-value", trans = "log10",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(name = "Gene\nCount", range = c(2, 8)) +
    labs(title = title,
         x = "Gene Ratio", 
         y = "GO Terms") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 10))
  
  return(p)
}

# Function to create network plot
create_go_network_plot <- function(go_results, fold_changes = NULL, title = "GO Network", 
                                  selected_terms = NULL, top_n = 15) {
  
  if (is.null(go_results) || nrow(go_results) == 0) {
    return(NULL)
  }
  
  # Filter to selected terms or top terms
  if (!is.null(selected_terms)) {
    # Filter by description matching
    go_filtered <- go_results
    go_df <- as.data.frame(go_filtered)
    matching_terms <- go_df$Description[go_df$Description %in% selected_terms]
    
    if (length(matching_terms) > 0) {
      go_filtered@result <- go_df[go_df$Description %in% matching_terms, ]
    } else {
      # If no exact matches, use pattern matching
      pattern <- paste(selected_terms, collapse = "|")
      matching_rows <- grepl(pattern, go_df$Description, ignore.case = TRUE)
      go_filtered@result <- go_df[matching_rows, ]
    }
  } else {
    go_filtered <- go_results
    if (nrow(go_filtered) > top_n) {
      go_filtered@result <- head(as.data.frame(go_filtered), top_n)
    }
  }
  
  if (nrow(go_filtered@result) == 0) {
    cat("No terms found for network plot\n")
    return(NULL)
  }
  
  # Create network plot
  tryCatch({
    if (!is.null(fold_changes)) {
      p <- cnetplot(go_filtered, 
                   foldChange = fold_changes,
                   showCategory = min(nrow(go_filtered@result), 10),
                   layout = "kk") +
        scale_color_gradient2(low = "blue", mid = "white", high = "red",
                             name = "Log2FC") +
        theme_void() +
        labs(title = title)
    } else {
      p <- cnetplot(go_filtered, 
                   showCategory = min(nrow(go_filtered@result), 10),
                   layout = "kk") +
        theme_void() +
        labs(title = title)
    }
    
    return(p)
    
  }, error = function(e) {
    cat("Network plot creation failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to create comprehensive Excel report
create_go_excel_report <- function(go_results_list, output_file) {
  
  wb <- createWorkbook()
  
  # Create summary sheet
  summary_data <- data.frame(
    Analysis = character(),
    Total_Terms = integer(),
    Significant_Terms = integer(),
    Top_Term = character(),
    Top_Pvalue = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (analysis_name in names(go_results_list)) {
    go_result <- go_results_list[[analysis_name]]
    
    if (!is.null(go_result) && nrow(go_result) > 0) {
      go_df <- as.data.frame(go_result)
      
      # Add to summary
      summary_data <- rbind(summary_data, data.frame(
        Analysis = analysis_name,
        Total_Terms = nrow(go_df),
        Significant_Terms = sum(go_df$p.adjust < P_VALUE_CUTOFF),
        Top_Term = go_df$Description[1],
        Top_Pvalue = go_df$p.adjust[1],
        stringsAsFactors = FALSE
      ))
      
      # Add individual sheet
      addWorksheet(wb, analysis_name, zoom = 200)
      writeData(wb, analysis_name, go_df)
    }
  }
  
  # Add summary sheet
  if (nrow(summary_data) > 0) {
    addWorksheet(wb, "Summary", zoom = 200)
    writeData(wb, "Summary", summary_data)
  }
  
  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Excel report saved:", output_file, "\n")
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

cat("Starting GO enrichment analysis...\n")
cat("Job ID:", JOB_ID, "\n")
cat("Species:", SPECIES, "\n")

# Set up organism database
org_db <- setup_organism_db(SPECIES)
cat("Using organism database:", class(org_db)[1], "\n")

# Load DEG data
if (!file.exists(INPUT_DEG_FILE)) {
  stop("Input DEG file not found: ", INPUT_DEG_FILE)
}

deg_data <- readRDS(INPUT_DEG_FILE)
cat("Loaded DEG data with", nrow(deg_data), "genes\n")

# Validate DEG data structure
required_cols <- c("gene")
if (!all(required_cols %in% colnames(deg_data))) {
  stop("DEG data must contain 'gene' column")
}

# ============================================================================
# Prepare Gene Sets for Analysis
# ============================================================================

cat("Preparing gene sets for analysis...\n")

# Create different gene sets based on criteria
gene_sets <- list()

# All genes
gene_sets[["all_genes"]] <- deg_data$gene

# Genes by log2FC threshold (if available)
if ("avg_log2FC" %in% colnames(deg_data)) {
  gene_sets[["high_log2fc"]] <- deg_data$gene[abs(deg_data$avg_log2FC) > LOG2FC_THRESHOLD]
  gene_sets[["upregulated"]] <- deg_data$gene[deg_data$avg_log2FC > LOG2FC_THRESHOLD]
  gene_sets[["downregulated"]] <- deg_data$gene[deg_data$avg_log2FC < -LOG2FC_THRESHOLD]
}

# Genes by significance (if available)
if ("p_val_adj" %in% colnames(deg_data)) {
  gene_sets[["significant"]] <- deg_data$gene[deg_data$p_val_adj < P_VALUE_CUTOFF]
}

# By cluster (if available)
if ("cluster" %in% colnames(deg_data)) {
  clusters <- unique(deg_data$cluster)
  for (cluster in head(clusters, 5)) {  # Limit to top 5 clusters
    gene_sets[[paste0("cluster_", cluster)]] <- deg_data$gene[deg_data$cluster == cluster]
  }
}

# Remove empty gene sets
gene_sets <- gene_sets[sapply(gene_sets, length) > 0]

cat("Created", length(gene_sets), "gene sets:\n")
for (set_name in names(gene_sets)) {
  cat("  -", set_name, ":", length(gene_sets[[set_name]]), "genes\n")
}

# ============================================================================
# Perform GO Enrichment Analysis
# ============================================================================

cat("Performing GO enrichment analysis for each gene set...\n")

go_results_list <- list()
fold_change_lists <- list()

for (set_name in names(gene_sets)) {
  cat("\nAnalyzing gene set:", set_name, "\n")
  
  genes <- gene_sets[[set_name]]
  
  if (length(genes) < 5) {
    cat("Skipping", set_name, "- too few genes\n")
    next
  }
  
  # Convert gene IDs
  gene_df <- convert_gene_ids(genes, org_db)
  
  if (nrow(gene_df) == 0) {
    cat("No valid Entrez IDs for", set_name, "\n")
    next
  }
  
  # Perform GO enrichment
  go_result <- perform_go_enrichment(gene_df$ENTREZID, org_db, ONTOLOGY)
  
  if (!is.null(go_result) && nrow(go_result) > 0) {
    go_results_list[[set_name]] <- go_result
    
    # Prepare fold change data if available
    if ("avg_log2FC" %in% colnames(deg_data)) {
      # Create named vector of fold changes
      fc_data <- deg_data[deg_data$gene %in% gene_df$SYMBOL, ]
      if (nrow(fc_data) > 0) {
        fold_changes <- setNames(fc_data$avg_log2FC, fc_data$gene)
        fold_change_lists[[set_name]] <- fold_changes
      }
    }
  }
}

cat("\nSuccessfully analyzed", length(go_results_list), "gene sets\n")

# ============================================================================
# Create Visualizations
# ============================================================================

if (length(go_results_list) > 0) {
  
  cat("Creating visualizations...\n")
  
  # Define output prefix
  output_prefix <- paste0("GO_", SPECIES, "_", JOB_ID)
  
  # Create plots for each analysis
  for (set_name in names(go_results_list)) {
    go_result <- go_results_list[[set_name]]
    fold_changes <- fold_change_lists[[set_name]]
    
    cat("Creating plots for:", set_name, "\n")
    
    # Bar plot
    p_bar <- create_go_barplot(go_result, 
                              title = paste("GO Enrichment -", set_name),
                              top_n = TOP_N_TERMS)
    
    if (!is.null(p_bar)) {
      bar_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_", set_name, "_barplot.pdf"))
      ggsave(bar_file, p_bar, width = 12, height = 8, dpi = 300)
    }
    
    # Dot plot
    p_dot <- create_go_dotplot(go_result, 
                              title = paste("GO Enrichment -", set_name),
                              top_n = TOP_N_TERMS)
    
    if (!is.null(p_dot)) {
      dot_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_", set_name, "_dotplot.pdf"))
      ggsave(dot_file, p_dot, width = 12, height = 8, dpi = 300)
    }
    
    # Network plot (if fold changes available)
    if (!is.null(fold_changes) && length(fold_changes) > 0) {
      p_net <- create_go_network_plot(go_result, fold_changes,
                                     title = paste("GO Network -", set_name))
      
      if (!is.null(p_net)) {
        net_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_", set_name, "_network.pdf"))
        ggsave(net_file, p_net, width = 12, height = 10, dpi = 300)
      }
    }
  }
  
  # Create combined summary plot
  if (length(go_results_list) > 1) {
    tryCatch({
      # Create comparison plot
      term_counts <- sapply(go_results_list, function(x) nrow(x))
      
      summary_df <- data.frame(
        Analysis = names(term_counts),
        Term_Count = term_counts
      )
      
      p_summary <- ggplot(summary_df, aes(x = reorder(Analysis, Term_Count), y = Term_Count)) +
        geom_col(fill = "steelblue", alpha = 0.7) +
        coord_flip() +
        labs(title = "GO Enrichment Summary",
             x = "Analysis", y = "Number of Significant Terms") +
        theme_minimal()
      
      summary_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_summary.pdf"))
      ggsave(summary_file, p_summary, width = 10, height = 6, dpi = 300)
      
    }, error = function(e) {
      cat("Error creating summary plot:", e$message, "\n")
    })
  }
}

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\n")

# Create Excel report
excel_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_results.xlsx"))
create_go_excel_report(go_results_list, excel_file)

# Save individual results as RDS
for (set_name in names(go_results_list)) {
  rds_file <- file.path(DATA_DIR, paste0(output_prefix, "_", set_name, ".rds"))
  saveRDS(go_results_list[[set_name]], rds_file)
}

# Save combined results
if (length(go_results_list) > 0) {
  combined_rds <- file.path(DATA_DIR, paste0(output_prefix, "_all_results.rds"))
  saveRDS(go_results_list, combined_rds)
}

# ============================================================================
# Generate Report
# ============================================================================

cat("Generating analysis report...\n")

n_gene_sets <- length(gene_sets)
n_successful <- length(go_results_list)

report_text <- paste0(
  "GO Enrichment Analysis Report\n",
  "============================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Job ID: ", JOB_ID, "\n",
  "Species: ", SPECIES, "\n",
  "Ontology: ", ONTOLOGY, "\n\n",
  "Parameters:\n",
  "  - P-value cutoff: ", P_VALUE_CUTOFF, "\n",
  "  - Q-value cutoff: ", Q_VALUE_CUTOFF, "\n",
  "  - Min gene set size: ", MIN_GS_SIZE, "\n",
  "  - Max gene set size: ", MAX_GS_SIZE, "\n",
  "  - Log2FC threshold: ", LOG2FC_THRESHOLD, "\n\n",
  "Input Data:\n",
  "  - Total genes: ", nrow(deg_data), "\n",
  "  - Gene sets created: ", n_gene_sets, "\n",
  "  - Successful analyses: ", n_successful, "\n\n"
)

if (n_successful > 0) {
  report_text <- paste0(report_text, "Analysis Results:\n")
  
  for (set_name in names(go_results_list)) {
    go_result <- go_results_list[[set_name]]
    n_terms <- nrow(go_result)
    top_term <- as.data.frame(go_result)$Description[1]
    
    report_text <- paste0(report_text,
      "  - ", set_name, ": ", n_terms, " significant terms\n",
      "    Top term: ", top_term, "\n"
    )
  }
}

report_text <- paste0(report_text,
  "\nOutput Files:\n",
  "  - Excel report: ", basename(excel_file), "\n",
  "  - Individual RDS files: ", paste0(output_prefix, "_[analysis].rds"), "\n",
  "  - Plots: ", paste0(output_prefix, "_[analysis]_[type].pdf"), "\n\n",
  "Recommendations:\n",
  "1. Review top enriched terms for biological relevance\n",
  "2. Validate key pathways with literature\n",
  "3. Consider pathway analysis (KEGG, Reactome)\n",
  "4. Examine gene-term networks for functional modules\n",
  "5. Compare enrichment across different conditions\n\n",
  "Usage Examples:\n",
  "Rscript go_enrichment_analysis.R <deg_file.rds> <species> <job_id>\n",
  "Rscript go_enrichment_analysis.R deg_results.rds human analysis_v1\n"
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
cat("GO enrichment analysis completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\nFinal Summary:\n")
cat("- Gene sets analyzed:", n_gene_sets, "\n")
cat("- Successful analyses:", n_successful, "\n")
cat("- Species:", SPECIES, "\n")
cat("- Ontology:", ONTOLOGY, "\n")

if (n_successful > 0) {
  total_terms <- sum(sapply(go_results_list, nrow))
  cat("- Total significant terms:", total_terms, "\n")
  
  best_analysis <- names(go_results_list)[which.max(sapply(go_results_list, nrow))]
  cat("- Most enriched analysis:", best_analysis, 
      "(", nrow(go_results_list[[best_analysis]]), "terms )\n")
}

cat("\nNext steps:\n")
cat("1. Review enrichment results and visualizations\n")
cat("2. Investigate top biological processes\n")
cat("3. Validate findings with functional experiments\n")
cat("4. Consider pathway-level analysis\n")

cat("\nAnalysis completed successfully!\n")