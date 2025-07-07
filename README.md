# 3D Spatial Transcriptome of Primate Early Organogenesis - Analysis Code

This repository contains the bioinformatics analysis code supporting the manuscript "Reconstruction of a three-dimensional spatial transcriptome of primate early organogenesis". This study represents the first comprehensive 3D spatial transcriptomic atlas of cynomolgus monkey embryos at Carnegie Stages 9 and 10, capturing key events of early organogenesis.

## Project Overview

### Scientific Context
Early organogenesis is a critical phase of embryogenesis that establishes the foundation for organ development. This study utilized Stereo-seq spatial transcriptomics to create a comprehensive 3D molecular atlas of primate embryos during this crucial developmental window.

### Key Findings
- **13 major tissue clusters** identified across CS9 and CS10 embryos
- **Tissue-specific subclustering** revealing specialized cell populations in:
  - Heart tube development (5 subclusters)
  - Gut tube regionalization (5 subclusters) 
  - Neural tissue organization (8 subclusters)
  - Somitogenesis processes (6 subclusters)
- **Spatially-resolved gene regulatory networks** via SCENIC analysis
- **Cell-cell communication networks** during organogenesis
- **Developmental trajectory analysis** capturing temporal transitions

### Methodological Innovations
- 3D reconstruction using Spateo alignment of consecutive tissue sections
- Integration of spatial transcriptomics with published scRNA-seq datasets
- Spatial cell communication analysis using CellChat
- Multi-modal trajectory analysis combining RNA velocity and pseudotime

## Repository Structure

```
├── clustering_annotation/          # Major cluster identification and annotation
├── integration_analysis/           # Spatial-scRNA integration and correlation
├── subclustering/                  # Tissue-specific subclustering analysis
├── gene_regulatory_networks/       # SCENIC transcription factor analysis
├── trajectory_analysis/            # Developmental trajectory inference
├── cell_communication/             # CellChat ligand-receptor analysis
├── comparative_analysis/           # CS9 vs CS10 temporal comparisons
└── functional_analysis/            # GO enrichment and pathway analysis
```

## System Requirements

### Software Dependencies
- **R version**: 4.3.1 or higher
- **Key R packages**:
  - Seurat (spatial and single-cell analysis)
  - Harmony (batch correction and integration)
  - CellChat (cell-cell communication)
  - SCENIC (gene regulatory networks)
  - Monocle3 (trajectory analysis)
  - velocyto.R (RNA velocity)
  - clusterProfiler (functional enrichment)

### Computing Environment
- **Memory**: Minimum 64GB RAM recommended for large spatial datasets
- **Storage**: At least 100GB free space for intermediate files
- **Environment**: Conda environment with R 4.3.1 recommended

## Data Requirements

### Input Data Formats
- **Spatial transcriptomics**: Stereo-seq bin50 resolution data
- **Reference datasets**: Published scRNA-seq datasets for integration
- **Metadata**: Cell type annotations and developmental stage information

## Usage Examples

### Basic Workflow
1. **Clustering and Annotation**
   ```r
   source("clustering_annotation/initial_clustering.R")
   ```

2. **Integration Analysis**
   ```r
   source("integration_analysis/scrna_spatial_integration.R")
   ```

3. **Tissue-Specific Analysis**
   ```r
   source("subclustering/heart_tube_subclustering.R")
   ```

### Key Parameters
- **Clustering resolution**: 0.5-0.8 for major clusters, 0.3-0.6 for subclusters
- **Integration method**: Harmony for batch correction
- **Trajectory analysis**: RNA velocity + Monocle3 pseudotime

## Analysis Overview

### Clustering and Annotation
Identification of 13 major tissue clusters using UMAP dimensionality reduction and marker gene-based annotation with spatial validation.

### Integration Analysis  
Cross-modal integration of spatial transcriptomics with published scRNA-seq datasets using Harmony batch correction and correlation analysis.

### Subclustering Analysis
Tissue-specific subclustering to identify specialized cell populations within major lineages, focusing on heart, gut, neural, and somite development.

### Gene Regulatory Networks
SCENIC analysis to identify transcription factor regulons and regulatory networks governing tissue specification and development.

### Trajectory Analysis
RNA velocity and pseudotime analysis to capture developmental transitions and cell fate decisions during organogenesis.

### Cell Communication
CellChat analysis to map ligand-receptor interactions and signaling pathways coordinating tissue development.

### Comparative Analysis
Temporal comparison between CS9 and CS10 to identify stage-specific gene expression changes and developmental progression.

### Functional Analysis
Gene Ontology enrichment and pathway analysis to understand functional significance of identified gene signatures.

## Citation

If you use this code in your research, please cite:

```
[Manuscript citation will be added upon publication]
```

## Data Availability

The spatial transcriptomics datasets supporting this study are available at:
- [Data repository information will be added upon publication]

## Contact

For questions about the analysis methods or code, please submit an [issue on GitHub](../../issues). This includes:
- Technical problems and bug reports
- Feature requests and suggestions
- Questions about implementation
- Scientific inquiries about the analysis methods

## License

This code is provided under the MIT License.

