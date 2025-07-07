#!/usr/bin/env Rscript

# ============================================================================
# Comprehensive Pathway Enrichment Analysis
# ============================================================================
# 
# This script performs comprehensive pathway enrichment analysis including
# KEGG pathways, Reactome pathways, and custom pathway databases for
# single-cell and spatial transcriptomics data.
#
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(DOSE)
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
  library(pathview)
})

# Parameters
JOB_ID <- "pathway_analysis"
SPECIES <- "human"                    # "human" or "mouse"
PATHWAY_DATABASES <- c("KEGG", "Reactome", "WikiPathways")  # Pathway databases to use
P_VALUE_CUTOFF <- 0.05               # P-value cutoff for significance
Q_VALUE_CUTOFF <- 0.05               # Q-value cutoff for significance
MIN_GS_SIZE <- 10                    # Minimum gene set size
MAX_GS_SIZE <- 500                   # Maximum gene set size
P_ADJUST_METHOD <- "fdr"             # P-value adjustment method
LOG2FC_THRESHOLD <- 0.5              # Log2FC threshold for gene filtering
TOP_N_PATHWAYS <- 20                 # Number of top pathways to show

# Custom pathway categories for developmental biology
DEVELOPMENTAL_PATHWAYS <- c(
  "embryonic", "development", "morphogenesis", "differentiation",
  "neural", "cardiac", "mesoderm", "endoderm", "ectoderm",
  "signaling", "transcription", "cell cycle", "apoptosis"
)

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

# Function to set up organism database and species code
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
    return(gene_df)
    
  }, error = function(e) {
    stop("Gene ID conversion failed: ", e$message)
  })
}

# Function to perform KEGG pathway enrichment
perform_kegg_enrichment <- function(entrez_ids, organism_code) {
  
  cat("Performing KEGG pathway enrichment...\n")
  
  if (length(entrez_ids) == 0) {
    return(NULL)
  }
  
  tryCatch({
    kegg_results <- enrichKEGG(
      gene = entrez_ids,
      organism = organism_code,
      keyType = "kegg",
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = P_ADJUST_METHOD,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      qvalueCutoff = Q_VALUE_CUTOFF,
      use_internal_data = FALSE
    )
    
    if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
      # Make readable
      kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      cat("Found", nrow(kegg_results), "significant KEGG pathways\n")
      return(kegg_results)
    } else {
      cat("No significant KEGG pathways found\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("KEGG enrichment failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to perform Reactome pathway enrichment
perform_reactome_enrichment <- function(entrez_ids, organism) {
  
  cat("Performing Reactome pathway enrichment...\n")
  
  if (length(entrez_ids) == 0) {
    return(NULL)
  }
  
  tryCatch({
    reactome_results <- enrichPathway(
      gene = entrez_ids,
      organism = organism,
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = P_ADJUST_METHOD,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE,
      qvalueCutoff = Q_VALUE_CUTOFF,
      readable = TRUE
    )
    
    if (!is.null(reactome_results) && nrow(reactome_results) > 0) {
      cat("Found", nrow(reactome_results), "significant Reactome pathways\n")
      return(reactome_results)
    } else {
      cat("No significant Reactome pathways found\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("Reactome enrichment failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to perform WikiPathways enrichment
perform_wikipathways_enrichment <- function(entrez_ids, organism_code) {
  
  cat("Performing WikiPathways enrichment...\n")
  
  if (length(entrez_ids) == 0) {
    return(NULL)
  }
  
  tryCatch({
    # Convert organism code for WikiPathways
    wp_organism <- switch(organism_code,
                         "hsa" = "Homo sapiens",
                         "mmu" = "Mus musculus")
    
    wp_results <- enrichWP(
      gene = entrez_ids,
      organism = wp_organism,
      pvalueCutoff = P_VALUE_CUTOFF,
      pAdjustMethod = P_ADJUST_METHOD,
      minGSSize = MIN_GS_SIZE,
      maxGSSize = MAX_GS_SIZE
    )
    
    if (!is.null(wp_results) && nrow(wp_results) > 0) {
      cat("Found", nrow(wp_results), "significant WikiPathways\n")
      return(wp_results)
    } else {
      cat("No significant WikiPathways found\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("WikiPathways enrichment failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to filter pathways by developmental relevance
filter_developmental_pathways <- function(pathway_results, keywords = DEVELOPMENTAL_PATHWAYS) {
  
  if (is.null(pathway_results) || nrow(pathway_results) == 0) {
    return(NULL)
  }
  
  pathway_df <- as.data.frame(pathway_results)
  
  # Create pattern to match developmental terms
  pattern <- paste(keywords, collapse = "|")
  
  # Filter by description
  dev_pathways <- pathway_df[grepl(pattern, pathway_df$Description, ignore.case = TRUE), ]
  
  cat("Filtered to", nrow(dev_pathways), "developmental pathways\n")
  
  return(dev_pathways)
}

# Function to create pathway dot plot
create_pathway_dotplot <- function(pathway_results, title = "Pathway Enrichment", top_n = TOP_N_PATHWAYS) {
  
  if (is.null(pathway_results) || nrow(pathway_results) == 0) {
    return(NULL)
  }
  
  pathway_df <- as.data.frame(pathway_results)
  
  # Select top pathways
  if (nrow(pathway_df) > top_n) {
    pathway_df <- head(pathway_df, top_n)
  }
  
  # Calculate gene ratio
  pathway_df$GeneRatio_numeric <- sapply(pathway_df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  
  # Create plot
  p <- pathway_df %>%
    mutate(Description = str_wrap(Description, width = 60)) %>%
    mutate(Description = fct_reorder(Description, Count)) %>%
    ggplot(aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.7) +
    scale_color_viridis_c(name = "Adjusted\nP-value", trans = "log10",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    labs(title = title,
         x = "Gene Ratio", 
         y = "Pathways") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 9))
  
  return(p)
}

# Function to create pathway bar plot
create_pathway_barplot <- function(pathway_results, title = "Pathway Enrichment", top_n = TOP_N_PATHWAYS) {
  
  if (is.null(pathway_results) || nrow(pathway_results) == 0) {
    return(NULL)
  }
  
  pathway_df <- as.data.frame(pathway_results)
  
  # Select top pathways
  if (nrow(pathway_df) > top_n) {
    pathway_df <- head(pathway_df, top_n)
  }
  
  # Create plot
  p <- pathway_df %>%
    mutate(Description = str_wrap(Description, width = 60)) %>%
    mutate(Description = fct_reorder(Description, Count)) %>%
    ggplot(aes(x = Description, y = Count)) +
    geom_col(aes(fill = p.adjust), alpha = 0.8) +
    scale_fill_viridis_c(name = "Adjusted\nP-value", trans = "log10",
                        guide = guide_colorbar(reverse = TRUE)) +
    coord_flip() +
    labs(title = title,
         x = "Pathways", 
         y = "Gene Count") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 9))
  
  return(p)
}

# Function to create pathway network plot
create_pathway_network <- function(pathway_results, fold_changes = NULL, title = "Pathway Network") {
  
  if (is.null(pathway_results) || nrow(pathway_results) == 0) {
    return(NULL)
  }
  
  # Limit to top pathways for visualization
  if (nrow(pathway_results) > 15) {
    pathway_results@result <- head(as.data.frame(pathway_results), 15)
  }
  
  tryCatch({
    if (!is.null(fold_changes) && length(fold_changes) > 0) {
      p <- cnetplot(pathway_results, 
                   foldChange = fold_changes,
                   showCategory = min(nrow(pathway_results), 10),
                   layout = "kk",
                   node_label = "category") +
        scale_color_gradient2(low = "blue", mid = "white", high = "red",
                             name = "Log2FC") +
        theme_void() +
        labs(title = title)
    } else {
      p <- cnetplot(pathway_results, 
                   showCategory = min(nrow(pathway_results), 10),
                   layout = "kk",
                   node_label = "category") +
        theme_void() +
        labs(title = title)
    }
    
    return(p)
    
  }, error = function(e) {
    cat("Network plot creation failed:", e$message, "\n")
    return(NULL)
  })
}

# Function to create comprehensive pathway comparison
create_pathway_comparison <- function(pathway_results_list) {
  
  if (length(pathway_results_list) == 0) {
    return(NULL)
  }
  
  # Extract pathway counts
  pathway_counts <- sapply(pathway_results_list, function(x) {
    if (is.null(x)) return(0)
    nrow(x)
  })
  
  # Create comparison data frame
  comparison_df <- data.frame(
    Database = names(pathway_counts),
    Pathways = pathway_counts
  )
  
  # Create plot
  p <- ggplot(comparison_df, aes(x = reorder(Database, Pathways), y = Pathways)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = Pathways), hjust = -0.1) +
    coord_flip() +
    labs(title = "Pathway Database Comparison",
         x = "Database", 
         y = "Number of Significant Pathways") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# Function to save KEGG pathway view
save_kegg_pathview <- function(pathway_id, fold_changes, output_dir, organism_code) {
  
  if (is.null(fold_changes) || length(fold_changes) == 0) {
    return(NULL)
  }
  
  tryCatch({
    # Convert gene symbols to entrez IDs for pathview
    gene_data <- fold_changes
    names(gene_data) <- names(fold_changes)
    
    # Create pathway view
    pathview(gene.data = gene_data,
             pathway.id = pathway_id,
             species = organism_code,
             out.suffix = paste0("_", JOB_ID),
             kegg.dir = output_dir)
    
    cat("Saved pathway view for", pathway_id, "\n")
    
  }, error = function(e) {
    cat("Pathview creation failed for", pathway_id, ":", e$message, "\n")
  })
}

# ============================================================================
# Main Analysis Pipeline
# ============================================================================

cat("Starting pathway enrichment analysis...\n")
cat("Job ID:", JOB_ID, "\n")
cat("Species:", SPECIES, "\n")

# Set up organism information
org_info <- setup_organism_info(SPECIES)
cat("Using organism:", org_info$kegg_organism, "\n")

# Load DEG data
if (!file.exists(INPUT_DEG_FILE)) {
  stop("Input DEG file not found: ", INPUT_DEG_FILE)
}

deg_data <- readRDS(INPUT_DEG_FILE)
cat("Loaded DEG data with", nrow(deg_data), "genes\n")

# Validate DEG data structure
if (!"gene" %in% colnames(deg_data)) {
  stop("DEG data must contain 'gene' column")
}

# ============================================================================
# Prepare Gene Sets
# ============================================================================

cat("Preparing gene sets...\n")

# Create different gene sets
gene_sets <- list()

# All genes
gene_sets[["all_genes"]] <- deg_data$gene

# Genes by log2FC (if available)
if ("avg_log2FC" %in% colnames(deg_data)) {
  gene_sets[["upregulated"]] <- deg_data$gene[deg_data$avg_log2FC > LOG2FC_THRESHOLD]
  gene_sets[["downregulated"]] <- deg_data$gene[deg_data$avg_log2FC < -LOG2FC_THRESHOLD]
  gene_sets[["high_fc"]] <- deg_data$gene[abs(deg_data$avg_log2FC) > LOG2FC_THRESHOLD]
}

# Significant genes (if available)
if ("p_val_adj" %in% colnames(deg_data)) {
  gene_sets[["significant"]] <- deg_data$gene[deg_data$p_val_adj < P_VALUE_CUTOFF]
}

# Remove empty sets
gene_sets <- gene_sets[sapply(gene_sets, length) > 5]

cat("Created", length(gene_sets), "gene sets\n")

# ============================================================================
# Perform Pathway Enrichment Analysis
# ============================================================================

cat("Performing pathway enrichment analysis...\n")

all_pathway_results <- list()

for (set_name in names(gene_sets)) {
  cat("\nAnalyzing gene set:", set_name, "\n")
  
  genes <- gene_sets[[set_name]]
  
  # Convert gene IDs
  gene_df <- convert_gene_ids(genes, org_info$org_db)
  
  if (nrow(gene_df) == 0) {
    cat("No valid Entrez IDs for", set_name, "\n")
    next
  }
  
  entrez_ids <- gene_df$ENTREZID
  
  # Perform enrichment for each database
  if ("KEGG" %in% PATHWAY_DATABASES) {
    kegg_result <- perform_kegg_enrichment(entrez_ids, org_info$kegg_organism)
    if (!is.null(kegg_result)) {
      all_pathway_results[[paste0(set_name, "_KEGG")]] <- kegg_result
    }
  }
  
  if ("Reactome" %in% PATHWAY_DATABASES) {
    reactome_result <- perform_reactome_enrichment(entrez_ids, org_info$reactome_organism)
    if (!is.null(reactome_result)) {
      all_pathway_results[[paste0(set_name, "_Reactome")]] <- reactome_result
    }
  }
  
  if ("WikiPathways" %in% PATHWAY_DATABASES) {
    wp_result <- perform_wikipathways_enrichment(entrez_ids, org_info$kegg_organism)
    if (!is.null(wp_result)) {
      all_pathway_results[[paste0(set_name, "_WikiPathways")]] <- wp_result
    }
  }
}

cat("\nSuccessfully analyzed", length(all_pathway_results), "pathway sets\n")

# ============================================================================
# Filter for Developmental Pathways
# ============================================================================

cat("Filtering for developmental pathways...\n")

dev_pathway_results <- list()

for (result_name in names(all_pathway_results)) {
  pathway_result <- all_pathway_results[[result_name]]
  dev_pathways <- filter_developmental_pathways(pathway_result)
  
  if (!is.null(dev_pathways) && nrow(dev_pathways) > 0) {
    dev_pathway_results[[paste0(result_name, "_dev")]] <- dev_pathways
  }
}

cat("Found", length(dev_pathway_results), "developmental pathway sets\n")

# ============================================================================
# Create Visualizations
# ============================================================================

if (length(all_pathway_results) > 0) {
  
  cat("Creating visualizations...\n")
  
  output_prefix <- paste0("pathway_", SPECIES, "_", JOB_ID)
  
  # Create plots for each analysis
  for (result_name in names(all_pathway_results)) {
    pathway_result <- all_pathway_results[[result_name]]
    
    cat("Creating plots for:", result_name, "\n")
    
    # Dot plot
    p_dot <- create_pathway_dotplot(pathway_result, 
                                   title = paste("Pathway Enrichment -", result_name))
    
    if (!is.null(p_dot)) {
      dot_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_", result_name, "_dotplot.pdf"))
      ggsave(dot_file, p_dot, width = 14, height = 10, dpi = 300)
    }
    
    # Bar plot
    p_bar <- create_pathway_barplot(pathway_result, 
                                   title = paste("Pathway Enrichment -", result_name))
    
    if (!is.null(p_bar)) {
      bar_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_", result_name, "_barplot.pdf"))
      ggsave(bar_file, p_bar, width = 14, height = 10, dpi = 300)
    }
    
    # Network plot (if fold changes available)
    if ("avg_log2FC" %in% colnames(deg_data)) {
      # Prepare fold changes
      fc_data <- deg_data[!is.na(deg_data$avg_log2FC), ]
      fold_changes <- setNames(fc_data$avg_log2FC, fc_data$gene)
      
      p_net <- create_pathway_network(pathway_result, fold_changes,
                                     title = paste("Pathway Network -", result_name))
      
      if (!is.null(p_net)) {
        net_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_", result_name, "_network.pdf"))
        ggsave(net_file, p_net, width = 12, height = 10, dpi = 300)
      }
    }
  }
  
  # Create database comparison plot
  # Group results by database
  db_results <- list()
  for (db in PATHWAY_DATABASES) {
    db_results[[db]] <- all_pathway_results[grepl(db, names(all_pathway_results))]
  }
  
  p_comparison <- create_pathway_comparison(lapply(db_results, function(x) {
    if (length(x) > 0) x[[1]] else NULL
  }))
  
  if (!is.null(p_comparison)) {
    comp_file <- file.path(PLOTS_DIR, paste0(output_prefix, "_database_comparison.pdf"))
    ggsave(comp_file, p_comparison, width = 10, height = 6, dpi = 300)
  }
}

# ============================================================================
# Save Results
# ============================================================================

cat("Saving results...\n")

# Create Excel workbook
wb <- createWorkbook()

# Add summary sheet
summary_data <- data.frame(
  Analysis = character(),
  Database = character(),
  Pathways = integer(),
  Top_Pathway = character(),
  Top_Pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (result_name in names(all_pathway_results)) {
  pathway_result <- all_pathway_results[[result_name]]
  pathway_df <- as.data.frame(pathway_result)
  
  # Extract database name
  db_name <- str_extract(result_name, "(KEGG|Reactome|WikiPathways)")
  
  # Add to summary
  summary_data <- rbind(summary_data, data.frame(
    Analysis = str_remove(result_name, "_(KEGG|Reactome|WikiPathways)"),
    Database = db_name,
    Pathways = nrow(pathway_df),
    Top_Pathway = pathway_df$Description[1],
    Top_Pvalue = pathway_df$p.adjust[1],
    stringsAsFactors = FALSE
  ))
  
  # Add sheet for each analysis
  sheet_name <- substr(result_name, 1, 31)  # Excel sheet name limit
  addWorksheet(wb, sheet_name, zoom = 200)
  writeData(wb, sheet_name, pathway_df)
}

# Add summary sheet
if (nrow(summary_data) > 0) {
  addWorksheet(wb, "Summary", zoom = 200)
  writeData(wb, "Summary", summary_data)
}

# Save Excel file
excel_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_results.xlsx"))
saveWorkbook(wb, excel_file, overwrite = TRUE)

# Save individual RDS files
for (result_name in names(all_pathway_results)) {
  rds_file <- file.path(DATA_DIR, paste0(output_prefix, "_", result_name, ".rds"))
  saveRDS(all_pathway_results[[result_name]], rds_file)
}

# Save combined results
combined_file <- file.path(DATA_DIR, paste0(output_prefix, "_all_results.rds"))
saveRDS(all_pathway_results, combined_file)

# ============================================================================
# Generate Report
# ============================================================================

cat("Generating analysis report...\n")

n_analyses <- length(all_pathway_results)
total_pathways <- sum(sapply(all_pathway_results, nrow))

report_text <- paste0(
  "Pathway Enrichment Analysis Report\n",
  "=================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Job ID: ", JOB_ID, "\n",
  "Species: ", SPECIES, "\n",
  "Databases: ", paste(PATHWAY_DATABASES, collapse = ", "), "\n\n",
  "Parameters:\n",
  "  - P-value cutoff: ", P_VALUE_CUTOFF, "\n",
  "  - Q-value cutoff: ", Q_VALUE_CUTOFF, "\n",
  "  - Min gene set size: ", MIN_GS_SIZE, "\n",
  "  - Max gene set size: ", MAX_GS_SIZE, "\n",
  "  - Log2FC threshold: ", LOG2FC_THRESHOLD, "\n\n",
  "Results Summary:\n",
  "  - Total analyses: ", n_analyses, "\n",
  "  - Total significant pathways: ", total_pathways, "\n",
  "  - Gene sets analyzed: ", length(gene_sets), "\n\n"
)

if (n_analyses > 0) {
  report_text <- paste0(report_text, "Analysis Details:\n")
  
  for (result_name in names(all_pathway_results)) {
    pathway_result <- all_pathway_results[[result_name]]
    n_pathways <- nrow(pathway_result)
    top_pathway <- as.data.frame(pathway_result)$Description[1]
    
    report_text <- paste0(report_text,
      "  - ", result_name, ": ", n_pathways, " pathways\n",
      "    Top: ", top_pathway, "\n"
    )
  }
}

report_text <- paste0(report_text,
  "\nOutput Files:\n",
  "  - Excel report: ", basename(excel_file), "\n",
  "  - Combined RDS: ", basename(combined_file), "\n",
  "  - Plots: ", paste0(output_prefix, "_[analysis]_[type].pdf"), "\n\n",
  "Recommendations:\n",
  "1. Focus on pathways with high gene counts and low p-values\n",
  "2. Validate developmental pathways with literature\n",
  "3. Consider pathway crosstalk and interactions\n",
  "4. Examine pathway networks for functional modules\n",
  "5. Compare pathways across different conditions\n\n",
  "Usage Examples:\n",
  "Rscript pathway_enrichment_analysis.R <deg_file.rds> <species> <job_id>\n",
  "Rscript pathway_enrichment_analysis.R deg_results.rds human analysis_v1\n"
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
cat("Pathway enrichment analysis completed successfully!\n")
cat("Results saved in:", OUTPUT_DIR, "\n")
cat("Plots saved in:", PLOTS_DIR, "\n")

cat("\nFinal Summary:\n")
cat("- Analyses performed:", n_analyses, "\n")
cat("- Total pathways found:", total_pathways, "\n")
cat("- Databases used:", paste(PATHWAY_DATABASES, collapse = ", "), "\n")

if (n_analyses > 0) {
  best_analysis <- names(all_pathway_results)[which.max(sapply(all_pathway_results, nrow))]
  cat("- Most enriched analysis:", best_analysis, 
      "(", nrow(all_pathway_results[[best_analysis]]), "pathways )\n")
}

cat("\nNext steps:\n")
cat("1. Review pathway enrichment results\n")
cat("2. Investigate pathway networks and interactions\n")
cat("3. Validate key pathways functionally\n")
cat("4. Compare with literature and databases\n")

cat("\nAnalysis completed successfully!\n")