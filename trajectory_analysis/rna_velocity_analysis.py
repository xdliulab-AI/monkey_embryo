#!/usr/bin/env python3

"""
RNA Velocity Analysis with scVelo
=================================

This script performs RNA velocity analysis using scVelo to infer directional
cell state transitions and developmental trajectories. It integrates spliced
and unspliced RNA counts to predict future cell states.

"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad
import warnings
from pathlib import Path

# Visualization imports
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

# Suppress warnings
warnings.filterwarnings('ignore')

# Default parameters
DEFAULT_MIN_SHARED_COUNTS = 30
DEFAULT_N_TOP_GENES = 2000
DEFAULT_N_PCS = 30
DEFAULT_N_NEIGHBORS = 30
DEFAULT_MODE = 'dynamical'  # 'stochastic', 'dynamical', or 'steady_state'
DEFAULT_MIN_LIKELIHOOD = 0.1
DEFAULT_MIN_CONFIDENCE = 0.75

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="RNA Velocity Analysis with scVelo",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    parser.add_argument(
        "--velocyto_file",
        type=str,
        required=True,
        help="Path to Velocyto loom file with spliced/unspliced counts"
    )
    
    parser.add_argument(
        "--expression_file",
        type=str,
        help="Path to 10X expression matrix directory (optional)"
    )
    
    parser.add_argument(
        "--metadata_file",
        type=str,
        help="Path to cell metadata CSV file with barcodes and cell types"
    )
    
    parser.add_argument(
        "--coordinates_file",
        type=str,
        help="Path to coordinates CSV file with UMAP/spatial coordinates"
    )
    
    # Analysis parameters
    parser.add_argument(
        "--min_shared_counts",
        type=int,
        default=DEFAULT_MIN_SHARED_COUNTS,
        help="Minimum shared counts for filtering"
    )
    
    parser.add_argument(
        "--n_top_genes",
        type=int,
        default=DEFAULT_N_TOP_GENES,
        help="Number of top variable genes"
    )
    
    parser.add_argument(
        "--n_pcs",
        type=int,
        default=DEFAULT_N_PCS,
        help="Number of principal components"
    )
    
    parser.add_argument(
        "--n_neighbors",
        type=int,
        default=DEFAULT_N_NEIGHBORS,
        help="Number of neighbors for graph construction"
    )
    
    parser.add_argument(
        "--mode",
        type=str,
        choices=['stochastic', 'dynamical', 'steady_state'],
        default=DEFAULT_MODE,
        help="Velocity estimation mode"
    )
    
    parser.add_argument(
        "--min_likelihood",
        type=float,
        default=DEFAULT_MIN_LIKELIHOOD,
        help="Minimum likelihood for dynamical mode"
    )
    
    parser.add_argument(
        "--min_confidence",
        type=float,
        default=DEFAULT_MIN_CONFIDENCE,
        help="Minimum confidence for filtering"
    )
    
    # Output options
    parser.add_argument(
        "--output_dir",
        type=str,
        default="../output",
        help="Output directory"
    )
    
    parser.add_argument(
        "--plots_dir",
        type=str,
        default="../plots",
        help="Plots directory"
    )
    
    parser.add_argument(
        "--sample_name",
        type=str,
        default="velocity_analysis",
        help="Sample name for output files"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=-1,
        help="Number of parallel jobs (-1 for all cores)"
    )
    
    return parser.parse_args()

def setup_scvelo():
    """Setup scVelo settings."""
    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.set_figure_params('scvelo')
    
    # Set scanpy settings
    sc.settings.verbosity = 3

def load_velocyto_data(velocyto_file, verbose=False):
    """Load Velocyto loom file."""
    if verbose:
        print(f"Loading Velocyto data from: {velocyto_file}")
    
    if not os.path.exists(velocyto_file):
        raise FileNotFoundError(f"Velocyto file not found: {velocyto_file}")
    
    # Load loom file
    adata = scv.read(velocyto_file, cache=True)
    
    # Make variable names unique
    adata.var_names_make_unique()
    
    if verbose:
        print(f"Loaded data: {adata.n_obs} cells x {adata.n_vars} genes")
        print(f"Layers available: {list(adata.layers.keys())}")
    
    return adata

def load_expression_data(expression_file, verbose=False):
    """Load 10X expression data (optional)."""
    if expression_file and os.path.exists(expression_file):
        if verbose:
            print(f"Loading expression data from: {expression_file}")
        
        expr_data = sc.read_10x_mtx(
            expression_file,
            var_names='gene_symbols',
            cache=True
        )
        expr_data.var_names_make_unique()
        
        if verbose:
            print(f"Expression data: {expr_data.n_obs} cells x {expr_data.n_vars} genes")
        
        return expr_data
    
    return None

def load_metadata(metadata_file, verbose=False):
    """Load cell metadata."""
    if metadata_file and os.path.exists(metadata_file):
        if verbose:
            print(f"Loading metadata from: {metadata_file}")
        
        metadata = pd.read_csv(metadata_file, index_col='Barcode')
        
        if verbose:
            print(f"Metadata loaded: {len(metadata)} cells")
            print(f"Columns: {list(metadata.columns)}")
        
        return metadata
    
    return None

def load_coordinates(coordinates_file, verbose=False):
    """Load coordinate data (UMAP or spatial)."""
    if coordinates_file and os.path.exists(coordinates_file):
        if verbose:
            print(f"Loading coordinates from: {coordinates_file}")
        
        coords = pd.read_csv(coordinates_file, index_col='Barcode')
        
        if verbose:
            print(f"Coordinates loaded: {len(coords)} cells")
            print(f"Columns: {list(coords.columns)}")
        
        return coords
    
    return None

def integrate_metadata(adata, metadata, coordinates, verbose=False):
    """Integrate metadata and coordinates into AnnData object."""
    
    # Find overlapping barcodes
    overlapping_barcodes = set(adata.obs_names)
    
    if metadata is not None:
        overlapping_barcodes = overlapping_barcodes.intersection(set(metadata.index))
    
    if coordinates is not None:
        overlapping_barcodes = overlapping_barcodes.intersection(set(coordinates.index))
    
    overlapping_barcodes = list(overlapping_barcodes)
    
    if verbose:
        print(f"Overlapping barcodes: {len(overlapping_barcodes)}")
    
    if len(overlapping_barcodes) == 0:
        raise ValueError("No overlapping barcodes found between datasets")
    
    # Subset to overlapping barcodes
    adata = adata[overlapping_barcodes].copy()
    
    # Add metadata
    if metadata is not None:
        for col in metadata.columns:
            adata.obs[col] = metadata.loc[adata.obs_names, col].values
        
        if verbose:
            print(f"Added metadata columns: {list(metadata.columns)}")
    
    # Add coordinates
    if coordinates is not None:
        coord_columns = [col for col in coordinates.columns 
                        if any(x in col.lower() for x in ['coordinate', 'umap', 'tsne', 'x', 'y'])]
        
        if len(coord_columns) >= 2:
            # Assume first two coordinate columns are X and Y
            coord_matrix = coordinates.loc[adata.obs_names, coord_columns[:2]].values
            adata.obsm['X_umap'] = coord_matrix
            
            if verbose:
                print(f"Added coordinates: {coord_columns[:2]}")
        else:
            if verbose:
                print("Warning: Could not identify coordinate columns")
    
    return adata

def preprocess_velocity_data(adata, min_shared_counts=30, n_top_genes=2000, verbose=False):
    """Preprocess velocity data."""
    if verbose:
        print("Preprocessing velocity data...")
        print(f"Before filtering: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # Filter and normalize
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=min_shared_counts,
        n_top_genes=n_top_genes
    )
    
    if verbose:
        print(f"After filtering: {adata.n_obs} cells x {adata.n_vars} genes")
    
    return adata

def compute_velocity(adata, mode='dynamical', n_pcs=30, n_neighbors=30, 
                    min_likelihood=0.1, n_jobs=-1, verbose=False):
    """Compute RNA velocity."""
    if verbose:
        print(f"Computing velocity using {mode} mode...")
    
    # Compute moments (first and second moments of spliced/unspliced)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    
    if mode == 'dynamical':
        # Recover dynamics
        scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
        
        # Filter genes based on likelihood
        scv.tl.velocity(adata, mode='dynamical')
        
        # Filter low-likelihood genes
        if min_likelihood > 0:
            gene_likelihood = adata.var['fit_likelihood']
            high_likelihood_genes = gene_likelihood > min_likelihood
            
            if verbose:
                print(f"Genes passing likelihood filter: {sum(high_likelihood_genes)}/{len(high_likelihood_genes)}")
            
            if sum(high_likelihood_genes) > 10:  # Ensure we have enough genes
                adata = adata[:, high_likelihood_genes].copy()
                scv.tl.velocity(adata, mode='dynamical')
    
    elif mode == 'stochastic':
        scv.tl.velocity(adata, mode='stochastic')
    
    elif mode == 'steady_state':
        scv.tl.velocity(adata, mode='steady_state')
    
    # Compute velocity graph
    scv.tl.velocity_graph(adata, n_jobs=n_jobs)
    
    if verbose:
        print("Velocity computation completed")
    
    return adata

def compute_additional_metrics(adata, verbose=False):
    """Compute additional velocity metrics."""
    if verbose:
        print("Computing additional velocity metrics...")
    
    # Compute velocity confidence
    scv.tl.velocity_confidence(adata)
    
    # Compute latent time (pseudotime)
    scv.tl.latent_time(adata)
    
    # Compute terminal states
    scv.tl.terminal_states(adata)
    
    # Compute driver genes
    try:
        scv.tl.rank_velocity_genes(adata, groupby='celltype')
    except:
        try:
            scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters')
        except:
            if verbose:
                print("Warning: Could not compute velocity genes (no grouping variable found)")
    
    if verbose:
        print("Additional metrics computed")
    
    return adata

def create_visualizations(adata, output_dir, sample_name, verbose=False):
    """Create velocity visualizations."""
    if verbose:
        print("Creating velocity visualizations...")
    
    plots_created = []
    
    # Determine grouping variable
    groupby_var = None
    for var in ['celltype', 'cell_type', 'seurat_clusters', 'cluster']:
        if var in adata.obs.columns:
            groupby_var = var
            break
    
    if groupby_var is None:
        groupby_var = 'louvain'  # fallback
        sc.tl.louvain(adata)
    
    # 1. Velocity embedding stream plot
    try:
        scv.pl.velocity_embedding_stream(
            adata,
            basis='umap',
            color=groupby_var,
            save=f'{sample_name}_velocity_stream.pdf',
            show=False,
            figsize=(8, 6)
        )
        plots_created.append("velocity_stream")
    except Exception as e:
        if verbose:
            print(f"Could not create stream plot: {e}")
    
    # 2. Velocity embedding (arrow plot)
    try:
        scv.pl.velocity_embedding(
            adata,
            basis='umap',
            color=groupby_var,
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            save=f'{sample_name}_velocity_arrows.pdf',
            show=False,
            figsize=(8, 6)
        )
        plots_created.append("velocity_arrows")
    except Exception as e:
        if verbose:
            print(f"Could not create arrow plot: {e}")
    
    # 3. Velocity confidence
    try:
        scv.pl.scatter(
            adata,
            color='velocity_confidence',
            cmap='viridis',
            save=f'{sample_name}_velocity_confidence.pdf',
            show=False,
            figsize=(8, 6)
        )
        plots_created.append("velocity_confidence")
    except Exception as e:
        if verbose:
            print(f"Could not create confidence plot: {e}")
    
    # 4. Latent time (velocity pseudotime)
    try:
        scv.pl.scatter(
            adata,
            color='latent_time',
            color_map='viridis',
            size=80,
            save=f'{sample_name}_latent_time.pdf',
            show=False,
            figsize=(8, 6)
        )
        plots_created.append("latent_time")
    except Exception as e:
        if verbose:
            print(f"Could not create latent time plot: {e}")
    
    # 5. Terminal states
    try:
        scv.pl.scatter(
            adata,
            color='terminal_states',
            save=f'{sample_name}_terminal_states.pdf',
            show=False,
            figsize=(8, 6)
        )
        plots_created.append("terminal_states")
    except Exception as e:
        if verbose:
            print(f"Could not create terminal states plot: {e}")
    
    # 6. Velocity gene ranking (if available)
    try:
        if 'rank_velocity_genes' in adata.uns.keys():
            scv.pl.velocity_genes(
                adata,
                save=f'{sample_name}_velocity_genes.pdf',
                show=False
            )
            plots_created.append("velocity_genes")
    except Exception as e:
        if verbose:
            print(f"Could not create velocity genes plot: {e}")
    
    # 7. Phase portraits for top velocity genes
    try:
        # Get top velocity genes
        if 'velocity_genes' in adata.var.columns:
            top_genes = adata.var.nlargest(4, 'velocity_genes').index.tolist()
        else:
            # Fallback: use highly variable genes
            top_genes = adata.var.nlargest(4, 'dispersions_norm').index.tolist()
        
        if len(top_genes) > 0:
            scv.pl.velocity(
                adata,
                top_genes,
                ncols=2,
                save=f'{sample_name}_phase_portraits.pdf',
                show=False,
                figsize=(12, 8)
            )
            plots_created.append("phase_portraits")
    except Exception as e:
        if verbose:
            print(f"Could not create phase portraits: {e}")
    
    if verbose:
        print(f"Created plots: {plots_created}")
    
    return plots_created

def save_results(adata, output_dir, sample_name, verbose=False):
    """Save velocity analysis results."""
    if verbose:
        print("Saving velocity analysis results...")
    
    # Save processed AnnData object
    adata_file = os.path.join(output_dir, f'{sample_name}_velocity_analysis.h5ad')
    adata.write(adata_file)
    
    if verbose:
        print(f"Saved AnnData object: {adata_file}")
    
    # Save velocity vectors
    if 'velocity' in adata.layers.keys():
        velocity_df = pd.DataFrame(
            adata.layers['velocity'],
            index=adata.obs_names,
            columns=adata.var_names
        )
        velocity_file = os.path.join(output_dir, f'{sample_name}_velocity_vectors.csv')
        velocity_df.to_csv(velocity_file)
        
        if verbose:
            print(f"Saved velocity vectors: {velocity_file}")
    
    # Save velocity metrics
    velocity_metrics = pd.DataFrame(index=adata.obs_names)
    
    for metric in ['velocity_confidence', 'latent_time', 'terminal_states']:
        if metric in adata.obs.columns:
            velocity_metrics[metric] = adata.obs[metric]
    
    if len(velocity_metrics.columns) > 0:
        metrics_file = os.path.join(output_dir, f'{sample_name}_velocity_metrics.csv')
        velocity_metrics.to_csv(metrics_file)
        
        if verbose:
            print(f"Saved velocity metrics: {metrics_file}")
    
    # Save gene-level results
    gene_results = adata.var.copy()
    gene_file = os.path.join(output_dir, f'{sample_name}_velocity_genes.csv')
    gene_results.to_csv(gene_file)
    
    if verbose:
        print(f"Saved gene results: {gene_file}")
    
    return {
        'adata': adata_file,
        'velocity_vectors': velocity_file if 'velocity' in adata.layers.keys() else None,
        'velocity_metrics': metrics_file if len(velocity_metrics.columns) > 0 else None,
        'gene_results': gene_file
    }

def generate_report(adata, args, output_files, plots_created, runtime):
    """Generate analysis report."""
    
    # Calculate summary statistics
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    
    # Velocity statistics
    velocity_stats = {}
    if 'velocity_confidence' in adata.obs.columns:
        velocity_stats['mean_confidence'] = adata.obs['velocity_confidence'].mean()
        velocity_stats['high_confidence_cells'] = (adata.obs['velocity_confidence'] > args.min_confidence).sum()
    
    if 'latent_time' in adata.obs.columns:
        velocity_stats['latent_time_range'] = (adata.obs['latent_time'].min(), adata.obs['latent_time'].max())
    
    # Cell type distribution
    cell_type_dist = {}
    for var in ['celltype', 'cell_type', 'seurat_clusters']:
        if var in adata.obs.columns:
            cell_type_dist = adata.obs[var].value_counts().to_dict()
            break
    
    # Create report
    report_lines = [
        "RNA Velocity Analysis Report",
        "===========================",
        "",
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Sample: {args.sample_name}",
        f"Mode: {args.mode}",
        "",
        "Input Data:",
        f"  - Velocyto file: {os.path.basename(args.velocyto_file)}",
        f"  - Cells analyzed: {n_cells:,}",
        f"  - Genes analyzed: {n_genes:,}",
        "",
        "Analysis Parameters:",
        f"  - Velocity mode: {args.mode}",
        f"  - Min shared counts: {args.min_shared_counts}",
        f"  - Top genes: {args.n_top_genes}",
        f"  - PCs: {args.n_pcs}",
        f"  - Neighbors: {args.n_neighbors}",
        f"  - Min likelihood: {args.min_likelihood}",
        "",
        "Results:",
        f"  - Runtime: {runtime:.1f} seconds",
    ]
    
    # Add velocity statistics
    if velocity_stats:
        report_lines.append("  - Velocity Statistics:")
        if 'mean_confidence' in velocity_stats:
            report_lines.append(f"    * Mean confidence: {velocity_stats['mean_confidence']:.3f}")
            report_lines.append(f"    * High confidence cells: {velocity_stats['high_confidence_cells']}")
        
        if 'latent_time_range' in velocity_stats:
            lt_min, lt_max = velocity_stats['latent_time_range']
            report_lines.append(f"    * Latent time range: [{lt_min:.3f}, {lt_max:.3f}]")
    
    # Add cell type distribution
    if cell_type_dist:
        report_lines.extend([
            "",
            "Cell Type Distribution:"
        ])
        for cell_type, count in list(cell_type_dist.items())[:10]:  # Top 10
            report_lines.append(f"  - {cell_type}: {count} cells")
    
    # Add output files
    report_lines.extend([
        "",
        "Output Files:"
    ])
    
    for file_type, file_path in output_files.items():
        if file_path:
            report_lines.append(f"  - {file_type}: {os.path.basename(file_path)}")
    
    # Add plots
    if plots_created:
        report_lines.extend([
            "",
            "Visualization Files:"
        ])
        for plot_type in plots_created:
            report_lines.append(f"  - {plot_type}: {args.sample_name}_{plot_type}.pdf")
    
    # Add recommendations
    report_lines.extend([
        "",
        "Recommendations:",
        "1. Review velocity confidence scores and filter low-confidence cells",
        "2. Validate velocity directions with known biological processes",
        "3. Analyze terminal states for endpoint identification",
        "4. Compare latent time with pseudotime from other methods",
        "5. Investigate velocity genes for regulatory mechanisms"
    ])
    
    # Save report
    report_path = os.path.join(args.output_dir, f'{args.sample_name}_velocity_report.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines))
    
    return report_path

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.plots_dir, exist_ok=True)
    
    # Setup scVelo
    setup_scvelo()
    
    if args.verbose:
        print("RNA Velocity Analysis with scVelo")
        print("=" * 50)
        print(f"Velocyto file: {args.velocyto_file}")
        print(f"Mode: {args.mode}")
        print(f"Sample: {args.sample_name}")
        print("=" * 50)
    
    try:
        import time
        start_time = time.time()
        
        # Load data
        adata = load_velocyto_data(args.velocyto_file, verbose=args.verbose)
        
        # Load additional data if provided
        expr_data = load_expression_data(args.expression_file, verbose=args.verbose)
        metadata = load_metadata(args.metadata_file, verbose=args.verbose)
        coordinates = load_coordinates(args.coordinates_file, verbose=args.verbose)
        
        # Integrate metadata
        if metadata is not None or coordinates is not None:
            adata = integrate_metadata(adata, metadata, coordinates, verbose=args.verbose)
        
        # Preprocess data
        adata = preprocess_velocity_data(
            adata,
            min_shared_counts=args.min_shared_counts,
            n_top_genes=args.n_top_genes,
            verbose=args.verbose
        )
        
        # Compute velocity
        adata = compute_velocity(
            adata,
            mode=args.mode,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            min_likelihood=args.min_likelihood,
            n_jobs=args.n_jobs,
            verbose=args.verbose
        )
        
        # Compute additional metrics
        adata = compute_additional_metrics(adata, verbose=args.verbose)
        
        # Create visualizations
        plots_created = create_visualizations(
            adata, args.plots_dir, args.sample_name, verbose=args.verbose
        )
        
        # Save results
        output_files = save_results(
            adata, args.output_dir, args.sample_name, verbose=args.verbose
        )
        
        # Calculate runtime
        runtime = time.time() - start_time
        
        # Generate report
        report_path = generate_report(adata, args, output_files, plots_created, runtime)
        
        if args.verbose:
            print(f"\nRNA velocity analysis completed successfully!")
            print(f"Runtime: {runtime:.1f} seconds")
            print(f"Results saved to: {args.output_dir}")
            print(f"Plots saved to: {args.plots_dir}")
            print(f"Report: {report_path}")
        
        return 0
        
    except Exception as e:
        print(f"Error during RNA velocity analysis: {str(e)}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())