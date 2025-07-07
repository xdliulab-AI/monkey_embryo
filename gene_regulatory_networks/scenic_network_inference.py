#!/usr/bin/env python3

"""
SCENIC Network Inference Script
==============================

This script performs gene regulatory network inference using GRNBoost2 or GENIE3
algorithms as part of the SCENIC pipeline for single-cell regulatory network analysis.

"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, genie3
from pyscenic.utils import modules_from_adjacencies
import loompy
import multiprocessing
import time
from pathlib import Path

# Default parameters
DEFAULT_N_WORKERS = 12
DEFAULT_ALGORITHM = "grnboost2"
DEFAULT_TF_FILE = "../data/allTFs_hg38.txt"
DEFAULT_INT_DIR = "../output/int"
DEFAULT_OUTPUT_FILE = "adj.tsv"

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="SCENIC Network Inference using GRNBoost2 or GENIE3",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--loom_file",
        type=str,
        required=True,
        help="Path to the loom file containing expression data"
    )
    
    parser.add_argument(
        "--tf_file",
        type=str,
        default=DEFAULT_TF_FILE,
        help="Path to transcription factor list file"
    )
    
    parser.add_argument(
        "--algorithm",
        type=str,
        choices=["grnboost2", "genie3"],
        default=DEFAULT_ALGORITHM,
        help="Network inference algorithm to use"
    )
    
    parser.add_argument(
        "--n_workers",
        type=int,
        default=DEFAULT_N_WORKERS,
        help="Number of workers for parallel processing"
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        default=DEFAULT_INT_DIR,
        help="Output directory for results"
    )
    
    parser.add_argument(
        "--output_file",
        type=str,
        default=DEFAULT_OUTPUT_FILE,
        help="Output filename for adjacency matrix"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed for reproducibility"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser.parse_args()

def validate_inputs(args):
    """Validate input files and parameters."""
    # Check loom file
    if not os.path.exists(args.loom_file):
        raise FileNotFoundError(f"Loom file not found: {args.loom_file}")
    
    # Check TF file
    if not os.path.exists(args.tf_file):
        raise FileNotFoundError(f"Transcription factor file not found: {args.tf_file}")
    
    # Check output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Validate number of workers
    max_workers = multiprocessing.cpu_count()
    if args.n_workers > max_workers:
        print(f"Warning: Requested {args.n_workers} workers, but only {max_workers} CPUs available")
        args.n_workers = max_workers
    
    return args

def load_expression_data(loom_file, verbose=False):
    """Load expression data from loom file."""
    if verbose:
        print(f"Loading expression data from: {loom_file}")
    
    # Load loom file
    with loompy.connect(loom_file, 'r') as ds:
        # Get gene names and cell barcodes
        gene_names = ds.ra.Gene if 'Gene' in ds.ra.keys() else [f"Gene_{i}" for i in range(ds.shape[0])]
        cell_names = ds.ca.CellID if 'CellID' in ds.ca.keys() else [f"Cell_{i}" for i in range(ds.shape[1])]
        
        # Load expression matrix (genes x cells)
        expression_matrix = ds[:, :]
    
    # Create DataFrame
    expr_df = pd.DataFrame(
        expression_matrix.T,  # Transpose to cells x genes
        index=cell_names,
        columns=gene_names
    )
    
    if verbose:
        print(f"Expression data shape: {expr_df.shape}")
        print(f"Genes: {len(gene_names)}")
        print(f"Cells: {len(cell_names)}")
    
    return expr_df

def load_transcription_factors(tf_file, expr_df, verbose=False):
    """Load transcription factor list and filter by available genes."""
    if verbose:
        print(f"Loading transcription factors from: {tf_file}")
    
    # Load TF names
    tf_names = load_tf_names(tf_file)
    
    # Filter TFs present in expression data
    available_genes = set(expr_df.columns)
    tf_names_filtered = [tf for tf in tf_names if tf in available_genes]
    
    if verbose:
        print(f"Total TFs in file: {len(tf_names)}")
        print(f"TFs present in data: {len(tf_names_filtered)}")
        print(f"TF coverage: {len(tf_names_filtered)/len(tf_names)*100:.1f}%")
    
    if len(tf_names_filtered) == 0:
        raise ValueError("No transcription factors found in expression data")
    
    return tf_names_filtered

def run_network_inference(expr_df, tf_names, algorithm="grnboost2", n_workers=12, seed=123, verbose=False):
    """Run network inference using specified algorithm."""
    if verbose:
        print(f"Running {algorithm.upper()} network inference...")
        print(f"Workers: {n_workers}")
        print(f"Transcription factors: {len(tf_names)}")
        print(f"Target genes: {expr_df.shape[1]}")
    
    # Set random seed
    np.random.seed(seed)
    
    # Record start time
    start_time = time.time()
    
    # Run network inference
    if algorithm.lower() == "grnboost2":
        adjacencies = grnboost2(
            expression_data=expr_df,
            tf_names=tf_names,
            verbose=verbose,
            seed=seed,
            client_or_address=n_workers
        )
    elif algorithm.lower() == "genie3":
        adjacencies = genie3(
            expression_data=expr_df,
            tf_names=tf_names,
            verbose=verbose,
            seed=seed,
            client_or_address=n_workers
        )
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    
    # Calculate runtime
    runtime = time.time() - start_time
    
    if verbose:
        print(f"Network inference completed in {runtime:.1f} seconds")
        print(f"Adjacencies generated: {len(adjacencies)}")
    
    return adjacencies

def save_adjacencies(adjacencies, output_path, verbose=False):
    """Save adjacency matrix to file."""
    if verbose:
        print(f"Saving adjacencies to: {output_path}")
    
    # Ensure proper column names
    if adjacencies.columns.tolist() != ['TF', 'target', 'importance']:
        adjacencies.columns = ['TF', 'target', 'importance']
    
    # Save to TSV file
    adjacencies.to_csv(output_path, sep='\t', index=False, header=True)
    
    if verbose:
        print(f"Adjacencies saved successfully")
        print(f"File size: {os.path.getsize(output_path) / 1024 / 1024:.1f} MB")

def generate_summary_report(args, expr_df, tf_names, adjacencies, runtime):
    """Generate a summary report of the network inference."""
    
    # Calculate statistics
    n_cells = expr_df.shape[0]
    n_genes = expr_df.shape[1]
    n_tfs = len(tf_names)
    n_adjacencies = len(adjacencies)
    mean_importance = adjacencies['importance'].mean()
    std_importance = adjacencies['importance'].std()
    
    # Create report
    report_lines = [
        "SCENIC Network Inference Report",
        "==============================",
        "",
        f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"Algorithm: {args.algorithm.upper()}",
        f"Random Seed: {args.seed}",
        "",
        "Input Data:",
        f"  - Expression file: {os.path.basename(args.loom_file)}",
        f"  - TF file: {os.path.basename(args.tf_file)}",
        f"  - Cells: {n_cells:,}",
        f"  - Genes: {n_genes:,}",
        f"  - Transcription factors: {n_tfs:,}",
        "",
        "Network Inference:",
        f"  - Algorithm: {args.algorithm.upper()}",
        f"  - Workers: {args.n_workers}",
        f"  - Runtime: {runtime:.1f} seconds",
        "",
        "Results:",
        f"  - Total adjacencies: {n_adjacencies:,}",
        f"  - Mean importance: {mean_importance:.6f}",
        f"  - Std importance: {std_importance:.6f}",
        f"  - Min importance: {adjacencies['importance'].min():.6f}",
        f"  - Max importance: {adjacencies['importance'].max():.6f}",
        "",
        "Top 10 Adjacencies by Importance:",
    ]
    
    # Add top adjacencies
    top_adjacencies = adjacencies.nlargest(10, 'importance')
    for _, row in top_adjacencies.iterrows():
        report_lines.append(f"  - {row['TF']} -> {row['target']}: {row['importance']:.6f}")
    
    report_lines.extend([
        "",
        "Output Files:",
        f"  - Adjacency matrix: {args.output_file}",
        f"  - Summary report: network_inference_report.txt",
        "",
        "Next Steps:",
        "1. Run SCENIC regulon analysis (scenic_regulon_analysis.R)",
        "2. Analyze regulon activity and cell states",
        "3. Perform downstream network analysis"
    ])
    
    # Save report
    report_path = os.path.join(args.output_dir, "network_inference_report.txt")
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    
    return report_path

def main():
    """Main function to run network inference."""
    # Parse arguments
    args = parse_arguments()
    
    # Validate inputs
    args = validate_inputs(args)
    
    # Print configuration
    if args.verbose:
        print("SCENIC Network Inference")
        print("=" * 50)
        print(f"Loom file: {args.loom_file}")
        print(f"TF file: {args.tf_file}")
        print(f"Algorithm: {args.algorithm}")
        print(f"Workers: {args.n_workers}")
        print(f"Output directory: {args.output_dir}")
        print(f"Seed: {args.seed}")
        print("=" * 50)
    
    try:
        # Load expression data
        expr_df = load_expression_data(args.loom_file, verbose=args.verbose)
        
        # Load transcription factors
        tf_names = load_transcription_factors(args.tf_file, expr_df, verbose=args.verbose)
        
        # Record start time
        start_time = time.time()
        
        # Run network inference
        adjacencies = run_network_inference(
            expr_df=expr_df,
            tf_names=tf_names,
            algorithm=args.algorithm,
            n_workers=args.n_workers,
            seed=args.seed,
            verbose=args.verbose
        )
        
        # Calculate total runtime
        total_runtime = time.time() - start_time
        
        # Save results
        output_path = os.path.join(args.output_dir, args.output_file)
        save_adjacencies(adjacencies, output_path, verbose=args.verbose)
        
        # Generate summary report
        report_path = generate_summary_report(args, expr_df, tf_names, adjacencies, total_runtime)
        
        # Print final summary
        if args.verbose:
            print("\nNetwork inference completed successfully!")
            print(f"Adjacencies: {len(adjacencies):,}")
            print(f"Runtime: {total_runtime:.1f} seconds")
            print(f"Results saved to: {output_path}")
            print(f"Report saved to: {report_path}")
        
        return 0
        
    except Exception as e:
        print(f"Error during network inference: {str(e)}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())