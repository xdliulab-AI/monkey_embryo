#!/usr/bin/env python3

# ============================================================================
# COMMOT Spatial Cell-Cell Communication Analysis
# ============================================================================
# 
# This script performs spatial cell-cell communication analysis using COMMOT
# to identify ligand-receptor interactions with spatial context in spatial
# transcriptomics data.
#
# ============================================================================

import os
import gc
import pickle
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import commot as ct

# Suppress warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# Parameters
SPECIES = "human"                    # "human" or "mouse" 
SIGNALING_TYPE = "Secreted Signaling"  # Signaling type for database
DISTANCE_THRESHOLD = 250             # Distance threshold for communication
MIN_CELL_PCT = 0.05                  # Minimum cell percentage for L-R filtering
HETEROMERIC = True                   # Include heteromeric complexes
PATHWAY_SUM = True                   # Sum pathway communications
N_NEIGHBORS = 30                     # Number of neighbors for spatial graph
SPATIAL_KEY = "spatial"              # Key for spatial coordinates

# Input/Output paths
INPUT_H5AD = "../data/annotated_adata_object.h5ad"
OUTPUT_DIR = "../output"
PLOTS_DIR = "../plots"

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="COMMOT spatial communication analysis")
    
    parser.add_argument("--input", "-i", type=str, default=INPUT_H5AD,
                       help="Input h5ad file path")
    parser.add_argument("--species", "-s", type=str, default=SPECIES,
                       choices=["human", "mouse"], help="Species")
    parser.add_argument("--signaling_type", "-t", type=str, default=SIGNALING_TYPE,
                       help="Signaling type for database")
    parser.add_argument("--distance", "-d", type=float, default=DISTANCE_THRESHOLD,
                       help="Distance threshold for communication")
    parser.add_argument("--min_pct", "-p", type=float, default=MIN_CELL_PCT,
                       help="Minimum cell percentage for L-R filtering")
    parser.add_argument("--output_dir", "-o", type=str, default=OUTPUT_DIR,
                       help="Output directory")
    
    return parser.parse_args()

def setup_directories(output_dir, plots_dir):
    """Create output directories if they don't exist."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    print(f"Plots directory: {plots_dir}")

def load_and_validate_data(input_file):
    """Load and validate input data."""
    print(f"Loading data from: {input_file}")
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Load AnnData object
    adata = sc.read_h5ad(input_file)
    print(f"Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes")
    
    # Validate required components
    cell_type_col = None
    for col in ["celltype", "preliminary_celltype", "leiden"]:
        if col in adata.obs.columns:
            cell_type_col = col
            break
    
    if cell_type_col is None:
        raise ValueError("No cell type annotation found in adata.obs")
    
    print(f"Using cell type column: {cell_type_col}")
    
    # Check for spatial coordinates
    spatial_coords = None
    coord_source = "none"
    
    # Try different coordinate sources
    if "spatial" in adata.obsm.keys():
        spatial_coords = adata.obsm["spatial"]
        coord_source = "spatial"
    elif "X_spatial" in adata.obsm.keys():
        spatial_coords = adata.obsm["X_spatial"]
        coord_source = "X_spatial"
    elif all(col in adata.obs.columns for col in ["x", "y"]):
        spatial_coords = adata.obs[["x", "y"]].values
        coord_source = "obs[x,y]"
    elif "X_umap" in adata.obsm.keys():
        spatial_coords = adata.obsm["X_umap"]
        coord_source = "X_umap"
        print("Warning: Using UMAP coordinates as spatial coordinates")
    
    if spatial_coords is None:
        raise ValueError("No spatial coordinates found")
    
    print(f"Using spatial coordinates from: {coord_source}")
    
    # Ensure spatial coordinates are in obsm
    if coord_source != "spatial":
        adata.obsm["spatial"] = spatial_coords
    
    # Set cell type as categorical
    adata.obs["celltype_analysis"] = adata.obs[cell_type_col].astype("category")
    
    print(f"Cell type distribution:")
    print(adata.obs["celltype_analysis"].value_counts())
    
    return adata, cell_type_col

def prepare_commot_database(species, signaling_type):
    """Prepare ligand-receptor database for COMMOT."""
    print(f"Loading {species} ligand-receptor database...")
    
    # Load the appropriate database
    if species.lower() == "human":
        df_lr = ct.pp.ligand_receptor_database(species='human', signaling_type=signaling_type)
    elif species.lower() == "mouse":
        df_lr = ct.pp.ligand_receptor_database(species='mouse', signaling_type=signaling_type)
    else:
        raise ValueError(f"Unsupported species: {species}")
    
    print(f"Loaded database with {len(df_lr)} ligand-receptor pairs")
    print(f"Signaling type: {signaling_type}")
    
    return df_lr

def filter_lr_database(df_lr, adata, min_cell_pct=0.05):
    """Filter ligand-receptor database based on expression."""
    print("Filtering ligand-receptor database...")
    
    # Filter database based on expression
    df_lr_filtered = ct.pp.filter_lr_database(df_lr, adata, min_cell_pct=min_cell_pct)
    
    print(f"Filtered to {len(df_lr_filtered)} ligand-receptor pairs")
    print(f"Minimum cell percentage threshold: {min_cell_pct}")
    
    return df_lr_filtered

def compute_spatial_communication(adata, df_lr_filtered, distance_threshold=250, 
                                heteromeric=True, pathway_sum=True):
    """Compute spatial communication using COMMOT."""
    print("Computing spatial communication...")
    
    # Run COMMOT spatial communication analysis
    ct.tl.spatial_communication(
        adata,
        database_name='cellchat',
        df_ligrec=df_lr_filtered,
        dis_thr=distance_threshold,
        heteromeric=heteromeric,
        pathway_sum=pathway_sum
    )
    
    print("Spatial communication analysis completed")
    
    # Check what was computed
    if 'commot-cellchat' in adata.obsp.keys():
        print("Communication scores computed and stored in adata.obsp['commot-cellchat']")
    
    if 'commot-cellchat-sum' in adata.obsm.keys():
        print("Pathway communication scores stored in adata.obsm['commot-cellchat-sum']")
    
    return adata

def analyze_communication_patterns(adata, cell_type_col="celltype_analysis"):
    """Analyze communication patterns from COMMOT results."""
    print("Analyzing communication patterns...")
    
    results = {}
    
    # Extract communication data
    if 'commot-cellchat' in adata.obsp.keys():
        comm_matrix = adata.obsp['commot-cellchat']
        results['communication_matrix'] = comm_matrix
        
        # Calculate summary statistics
        total_communication = np.sum(comm_matrix)
        mean_communication = np.mean(comm_matrix)
        
        results['total_communication'] = total_communication
        results['mean_communication'] = mean_communication
        
        print(f"Total communication strength: {total_communication:.2f}")
        print(f"Mean communication per cell pair: {mean_communication:.4f}")
    
    # Analyze by cell type
    if cell_type_col in adata.obs.columns:
        cell_types = adata.obs[cell_type_col].cat.categories
        n_types = len(cell_types)
        
        # Create cell type communication matrix
        comm_by_type = np.zeros((n_types, n_types))
        
        for i, ct1 in enumerate(cell_types):
            for j, ct2 in enumerate(cell_types):
                mask1 = adata.obs[cell_type_col] == ct1
                mask2 = adata.obs[cell_type_col] == ct2
                
                if 'commot-cellchat' in adata.obsp.keys():
                    submatrix = adata.obsp['commot-cellchat'][mask1, :][:, mask2]
                    comm_by_type[i, j] = np.mean(submatrix)
        
        # Create DataFrame for easier handling
        comm_df = pd.DataFrame(comm_by_type, 
                              index=cell_types, 
                              columns=cell_types)
        
        results['communication_by_celltype'] = comm_df
        
        # Find top communicating pairs
        comm_pairs = []
        for i, ct1 in enumerate(cell_types):
            for j, ct2 in enumerate(cell_types):
                if i != j:  # Exclude self-communication
                    comm_pairs.append({
                        'source': ct1,
                        'target': ct2,
                        'communication_strength': comm_by_type[i, j]
                    })
        
        comm_pairs_df = pd.DataFrame(comm_pairs)
        comm_pairs_df = comm_pairs_df.sort_values('communication_strength', ascending=False)
        
        results['top_communication_pairs'] = comm_pairs_df.head(20)
        
        print(f"Top 5 communicating cell type pairs:")
        for _, row in comm_pairs_df.head(5).iterrows():
            print(f"  {row['source']} -> {row['target']}: {row['communication_strength']:.4f}")
    
    return results

def create_visualizations(adata, results, output_prefix, plots_dir):
    """Create comprehensive visualizations."""
    print("Creating visualizations...")
    
    plots_created = []
    
    # 1. Communication heatmap by cell type
    if 'communication_by_celltype' in results:
        plt.figure(figsize=(10, 8))
        sns.heatmap(results['communication_by_celltype'], 
                   annot=True, cmap='viridis', fmt='.3f')
        plt.title('Cell-Cell Communication Strength by Cell Type')
        plt.xlabel('Target Cell Type')
        plt.ylabel('Source Cell Type')
        plt.tight_layout()
        
        plot_file = os.path.join(plots_dir, f"{output_prefix}_communication_heatmap.pdf")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        plots_created.append("communication_heatmap")
        print(f"Saved communication heatmap: {plot_file}")
    
    # 2. Top communication pairs bar plot
    if 'top_communication_pairs' in results:
        plt.figure(figsize=(12, 8))
        top_pairs = results['top_communication_pairs'].head(15)
        pair_labels = [f"{row['source']} → {row['target']}" 
                      for _, row in top_pairs.iterrows()]
        
        plt.barh(range(len(top_pairs)), top_pairs['communication_strength'])
        plt.yticks(range(len(top_pairs)), pair_labels)
        plt.xlabel('Communication Strength')
        plt.title('Top Cell Type Communication Pairs')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        
        plot_file = os.path.join(plots_dir, f"{output_prefix}_top_pairs.pdf")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        plots_created.append("top_pairs")
        print(f"Saved top pairs plot: {plot_file}")
    
    # 3. Spatial visualization if possible
    if 'spatial' in adata.obsm.keys():
        try:
            # Create spatial plot colored by total outgoing communication
            if 'commot-cellchat' in adata.obsp.keys():
                outgoing_comm = np.array(adata.obsp['commot-cellchat'].sum(axis=1)).flatten()
                adata.obs['outgoing_communication'] = outgoing_comm
                
                plt.figure(figsize=(10, 8))
                scatter = plt.scatter(adata.obsm['spatial'][:, 0], 
                                    adata.obsm['spatial'][:, 1],
                                    c=outgoing_comm, cmap='viridis', s=1)
                plt.colorbar(scatter, label='Outgoing Communication')
                plt.title('Spatial Distribution of Outgoing Communication')
                plt.xlabel('Spatial X')
                plt.ylabel('Spatial Y')
                plt.axis('equal')
                plt.tight_layout()
                
                plot_file = os.path.join(plots_dir, f"{output_prefix}_spatial_communication.pdf")
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plt.close()
                plots_created.append("spatial_communication")
                print(f"Saved spatial communication plot: {plot_file}")
        except Exception as e:
            print(f"Warning: Could not create spatial plot: {e}")
    
    return plots_created

def save_results(adata, results, output_prefix, output_dir):
    """Save analysis results."""
    print("Saving results...")
    
    files_saved = []
    
    # Save the processed AnnData object
    output_h5ad = os.path.join(output_dir, f"{output_prefix}_commot_results.h5ad")
    adata.write_h5ad(output_h5ad, compression='gzip')
    files_saved.append(output_h5ad)
    print(f"Saved processed data: {output_h5ad}")
    
    # Save communication results as CSV files
    if 'communication_by_celltype' in results:
        comm_file = os.path.join(output_dir, f"{output_prefix}_communication_by_celltype.csv")
        results['communication_by_celltype'].to_csv(comm_file)
        files_saved.append(comm_file)
        print(f"Saved communication matrix: {comm_file}")
    
    if 'top_communication_pairs' in results:
        pairs_file = os.path.join(output_dir, f"{output_prefix}_top_communication_pairs.csv")
        results['top_communication_pairs'].to_csv(pairs_file, index=False)
        files_saved.append(pairs_file)
        print(f"Saved top pairs: {pairs_file}")
    
    # Save summary statistics
    summary_stats = {
        'total_cells': adata.n_obs,
        'total_genes': adata.n_vars,
        'cell_types': len(adata.obs['celltype_analysis'].cat.categories),
        'total_communication': results.get('total_communication', 0),
        'mean_communication': results.get('mean_communication', 0)
    }
    
    summary_file = os.path.join(output_dir, f"{output_prefix}_summary_stats.txt")
    with open(summary_file, 'w') as f:
        for key, value in summary_stats.items():
            f.write(f"{key}: {value}\n")
    files_saved.append(summary_file)
    print(f"Saved summary statistics: {summary_file}")
    
    return files_saved

def generate_report(adata, results, output_prefix, output_dir, plots_created, files_saved):
    """Generate comprehensive analysis report."""
    print("Generating analysis report...")
    
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    n_cell_types = len(adata.obs['celltype_analysis'].cat.categories)
    
    report_text = f"""COMMOT Spatial Cell-Cell Communication Analysis Report
=====================================================

Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Species: {SPECIES}
Signaling Type: {SIGNALING_TYPE}

Data Summary:
  - Total cells: {n_cells:,}
  - Total genes: {n_genes:,}  
  - Cell types: {n_cell_types}

Analysis Parameters:
  - Distance threshold: {DISTANCE_THRESHOLD}
  - Minimum cell percentage: {MIN_CELL_PCT}
  - Heteromeric complexes: {HETEROMERIC}
  - Pathway summation: {PATHWAY_SUM}

Results Summary:"""
    
    if 'total_communication' in results:
        report_text += f"\n  - Total communication strength: {results['total_communication']:.2f}"
    if 'mean_communication' in results:
        report_text += f"\n  - Mean communication per pair: {results['mean_communication']:.4f}"
    
    # Add top communication pairs
    if 'top_communication_pairs' in results:
        report_text += f"\n\nTop 10 Communication Pairs:\n"
        top_pairs = results['top_communication_pairs'].head(10)
        for i, (_, row) in enumerate(top_pairs.iterrows(), 1):
            report_text += f"  {i:2d}. {row['source']} → {row['target']}: {row['communication_strength']:.4f}\n"
    
    # Add cell type information
    if 'celltype_analysis' in adata.obs.columns:
        report_text += f"\nCell Type Distribution:\n"
        cell_counts = adata.obs['celltype_analysis'].value_counts()
        for cell_type, count in cell_counts.items():
            pct = (count / n_cells) * 100
            report_text += f"  - {cell_type}: {count:,} cells ({pct:.1f}%)\n"
    
    # Add output files
    report_text += f"\nOutput Files:\n"
    for file_path in files_saved:
        filename = os.path.basename(file_path)
        report_text += f"  - {filename}\n"
    
    # Add visualization files
    if plots_created:
        report_text += f"\nVisualization Files:\n"
        for plot_type in plots_created:
            report_text += f"  - {output_prefix}_{plot_type}.pdf\n"
    
    report_text += f"""
Recommendations:
1. Examine spatial distribution of high-communication cell types
2. Validate key ligand-receptor pairs with expression analysis
3. Compare communication patterns across developmental stages
4. Investigate pathway-specific communication networks
5. Perform functional enrichment of communication pathways

Next Steps:
1. Analyze pathway-specific communication patterns
2. Investigate temporal changes in communication
3. Validate findings with experimental approaches
4. Compare with other communication analysis methods

Usage Examples:
python commot_spatial_communication.py --input data.h5ad --species human
python commot_spatial_communication.py --input data.h5ad --species mouse --distance 200
python commot_spatial_communication.py --input data.h5ad --signaling_type "ECM-Receptor"
"""
    
    report_file = os.path.join(output_dir, f"{output_prefix}_analysis_report.txt")
    with open(report_file, 'w') as f:
        f.write(report_text)
    
    print(f"Analysis report saved: {report_file}")
    return report_file

def main():
    """Main analysis pipeline."""
    print("COMMOT Spatial Cell-Cell Communication Analysis")
    print("=" * 50)
    
    # Parse arguments
    args = parse_arguments()
    
    # Setup directories
    setup_directories(args.output_dir, PLOTS_DIR)
    
    # Load and validate data
    adata, cell_type_col = load_and_validate_data(args.input)
    
    # Prepare database
    df_lr = prepare_commot_database(args.species, args.signaling_type)
    
    # Filter database
    df_lr_filtered = filter_lr_database(df_lr, adata, args.min_pct)
    
    # Compute spatial communication
    adata = compute_spatial_communication(
        adata, df_lr_filtered, 
        distance_threshold=args.distance,
        heteromeric=HETEROMERIC,
        pathway_sum=PATHWAY_SUM
    )
    
    # Analyze communication patterns
    results = analyze_communication_patterns(adata, cell_type_col="celltype_analysis")
    
    # Create output prefix
    output_prefix = f"commot_{args.species}_{args.signaling_type.replace(' ', '_')}"
    
    # Create visualizations
    plots_created = create_visualizations(adata, results, output_prefix, PLOTS_DIR)
    
    # Save results
    files_saved = save_results(adata, results, output_prefix, args.output_dir)
    
    # Generate report
    report_file = generate_report(adata, results, output_prefix, args.output_dir, 
                                plots_created, files_saved)
    
    # Final summary
    print("\n" + "=" * 50)
    print("COMMOT analysis completed successfully!")
    print(f"Results saved in: {args.output_dir}")
    print(f"Plots saved in: {PLOTS_DIR}")
    print(f"Report: {os.path.basename(report_file)}")
    
    print(f"\nSummary:")
    print(f"- Cells analyzed: {adata.n_obs:,}")
    print(f"- Cell types: {len(adata.obs['celltype_analysis'].cat.categories)}")
    print(f"- L-R pairs analyzed: {len(df_lr_filtered)}")
    
    if 'total_communication' in results:
        print(f"- Total communication: {results['total_communication']:.2f}")
    
    # Clean up
    gc.collect()

if __name__ == "__main__":
    main()