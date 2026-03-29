# Single-Cell RNA-seq Analysis Workflow

A comprehensive, reproducible workflow for analyzing single-cell RNA-seq data from raw counts to cell type annotation and trajectory analysis. Demonstrates both R (Seurat) and Python (Scanpy) approaches for maximum flexibility.

## Overview

This project showcases a complete single-cell RNA-seq analysis pipeline:

1. **Quality Control** — Cell filtering, gene filtering, doublet detection
2. **Normalization & Scaling** — Library size normalization, log transformation, scaling
3. **Dimensionality Reduction** — PCA, UMAP, t-SNE
4. **Clustering** — Leiden/Louvain clustering, cluster annotation
5. **Differential Expression** — Cluster-specific DE analysis
6. **Cell Type Annotation** — Marker-based and reference-based annotation
7. **Trajectory Analysis** — Pseudotime ordering, developmental trajectories
8. **Integration** — Multi-sample/multi-batch integration

## Skills Demonstrated

✅ **Single-Cell Tools:** Seurat (R), Scanpy (Python), CellRanger  
✅ **R Analysis:** tidyverse, ggplot2, Seurat, SingleCellExperiment  
✅ **Python Analysis:** pandas, NumPy, scikit-learn, Scanpy, matplotlib, seaborn  
✅ **Jupyter Notebooks:** Interactive exploratory analysis  
✅ **Data Integration:** Harmony, Seurat integration, batch correction  
✅ **Trajectory Analysis:** Monocle3, velocyto RNA velocity  
✅ **Visualization:** UMAP plots, violin plots, feature plots, heatmaps  
✅ **Reproducibility:** Conda environments, Docker containerization  
✅ **Version Control:** Git workflow with meaningful commits  

## Quick Start

### Prerequisites

- Docker (recommended) or Conda
- Git
- 16+ GB RAM for analysis
- 50+ GB disk space

### Installation

```bash
# Clone repository
git clone https://github.com/adamhoffman2155-hue/project-2-scrnaseq-analysis.git
cd project-2-scrnaseq-analysis

# Option 1: Using Docker
docker build -t scrnaseq-analysis .
docker run -it -v $(pwd):/workspace -p 8888:8888 scrnaseq-analysis bash

# Option 2: Using Conda
conda env create -f environment.yml
conda activate scrnaseq-analysis

# Start Jupyter
jupyter lab
```

## Project Structure

```
project-2-scrnaseq-analysis/
├── README.md                              # This file
├── Dockerfile                             # Container specification
├── environment.yml                        # Conda dependencies
├── notebooks/
│   ├── 01_qc_filtering.ipynb             # QC and filtering
│   ├── 02_normalization_scaling.ipynb    # Normalization workflow
│   ├── 03_dimensionality_reduction.ipynb # PCA, UMAP, t-SNE
│   ├── 04_clustering.ipynb               # Leiden clustering
│   ├── 05_cell_type_annotation.ipynb     # Marker-based annotation
│   ├── 06_differential_expression.ipynb  # Cluster DE analysis
│   ├── 07_trajectory_analysis.ipynb      # Pseudotime & velocity
│   └── 08_integration.ipynb              # Multi-sample integration
├── scripts/
│   ├── qc_metrics.py                     # QC calculation functions
│   ├── clustering.R                      # Seurat clustering wrapper
│   ├── annotation.py                     # Cell type annotation
│   ├── trajectory.R                      # Monocle3 trajectory analysis
│   └── visualization_utils.py            # Plotting utilities
├── data/
│   ├── raw/                              # Raw count matrices
│   ├── processed/                        # Processed data objects
│   └── metadata/                         # Cell metadata, annotations
├── results/
│   ├── qc/                               # QC reports
│   ├── figures/                          # Publication plots
│   ├── tables/                           # DE results, annotations
│   └── objects/                          # Seurat/AnnData objects
└── config/
    └── analysis_config.yaml              # Analysis parameters
```

## Workflow Details

### Phase 1: Quality Control

```
Raw Count Matrix (cells × genes)
    ↓
[Cell Filtering] → Remove low-quality cells (nUMI, nGenes, %MT)
    ↓
[Gene Filtering] → Remove lowly-expressed genes
    ↓
[Doublet Detection] → Remove doublets (DoubletFinder/Scrublet)
    ↓
Filtered Count Matrix
```

**Key Metrics:**
- nUMI (total UMI counts per cell)
- nGenes (number of detected genes)
- %MT (mitochondrial gene percentage)
- %RB (ribosomal gene percentage)

### Phase 2: Normalization & Scaling

```
Filtered Counts
    ↓
[Library Normalization] → Normalize to median library size
    ↓
[Log Transformation] → log(counts + 1)
    ↓
[Feature Selection] → Identify highly variable genes (HVGs)
    ↓
[Scaling] → Center and scale to unit variance
    ↓
Normalized & Scaled Matrix
```

### Phase 3: Dimensionality Reduction

```
Scaled Matrix
    ↓
[PCA] → Linear dimensionality reduction (50 PCs)
    ↓
[UMAP/t-SNE] → Non-linear visualization (2D)
    ↓
2D Embeddings for Visualization
```

### Phase 4: Clustering

```
PCA Embeddings
    ↓
[KNN Graph] → Build k-nearest neighbor graph
    ↓
[Leiden Algorithm] → Community detection clustering
    ↓
[Cluster Assignment] → Assign cells to clusters
    ↓
Cluster Labels
```

### Phase 5: Cell Type Annotation

```
Cluster-Specific Gene Expression
    ↓
[Marker Detection] → Find cluster-specific markers
    ↓
[Database Matching] → Compare to known markers (ImmGen, PanglaoDB)
    ↓
[Reference Mapping] → Map to reference single-cell datasets
    ↓
Cell Type Annotations
```

## Usage Examples

### Example 1: Run complete analysis (R/Seurat)

```bash
# Open Jupyter and run notebooks in order
jupyter lab

# Or run from command line
Rscript scripts/clustering.R --input data/raw/counts.mtx --output results/
```

### Example 2: Run Python analysis (Scanpy)

```python
import scanpy as sc
import pandas as pd

# Load data
adata = sc.read_h5ad('data/raw/counts.h5ad')

# QC
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 5000]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# PCA
sc.tl.pca(adata)

# UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Clustering
sc.tl.leiden(adata, resolution=0.5)

# Visualization
sc.pl.umap(adata, color='leiden', save='umap_clusters.pdf')
```

### Example 3: Cell type annotation

```python
# Load marker genes
markers = pd.read_csv('data/metadata/marker_genes.csv')

# Compute marker scores
for cell_type, genes in markers.groupby('cell_type'):
    sc.tl.score_genes(adata, genes['gene'].tolist(), score_name=cell_type)

# Assign cell types based on scores
adata.obs['cell_type'] = adata.obs[marker_cols].idxmax(axis=1)
```

## Input Data Formats

### Supported Formats

- **H5AD** (HDF5 Anndata) — recommended for Python
- **H5** (HDF5 10x Genomics) — standard CellRanger output
- **MTX** (Matrix Market) — sparse matrix format
- **CSV/TSV** — dense count matrices
- **Loom** — efficient single-cell format

### Example Metadata

```csv
cell_id,sample,condition,cell_type,n_genes,n_counts,pct_mt
cell_001,sample_1,control,T_cell,2500,15000,3.2
cell_002,sample_1,control,B_cell,3000,18000,2.8
cell_003,sample_1,treated,T_cell,2200,12000,4.1
```

## Output Files

### Main Results

- `results/objects/seurat_object.rds` — Complete Seurat object
- `results/objects/adata.h5ad` — Complete AnnData object
- `results/tables/cluster_markers.csv` — Cluster-specific markers
- `results/tables/cell_type_annotations.csv` — Cell type assignments
- `results/tables/de_results.csv` — Differential expression results

### Visualizations

- `results/figures/umap_clusters.pdf` — UMAP colored by cluster
- `results/figures/umap_celltypes.pdf` — UMAP colored by cell type
- `results/figures/marker_heatmap.pdf` — Top markers per cluster
- `results/figures/violin_plots.pdf` — Gene expression by cluster
- `results/figures/trajectory_plot.pdf` — Pseudotime trajectory

## Interpretation Guide

### UMAP Plot

- **Proximity:** Cells close in UMAP space are transcriptomically similar
- **Clusters:** Distinct regions represent cell populations
- **Colors:** Can represent clusters, cell types, or gene expression

### Marker Heatmap

- **Rows:** Top marker genes per cluster
- **Columns:** Individual cells (grouped by cluster)
- **Color intensity:** log-normalized expression level

### Violin Plots

- **X-axis:** Clusters or cell types
- **Y-axis:** Gene expression level
- **Shape:** Distribution of expression within cluster

## Performance Benchmarks

**Analysis Time** (10,000 cells, 20,000 genes, 4 cores):
- QC & Filtering: 2-5 minutes
- Normalization & Scaling: 5-10 minutes
- PCA: 2-5 minutes
- UMAP: 5-10 minutes
- Clustering: 2-5 minutes
- DE Analysis: 10-20 minutes
- **Total:** ~30-60 minutes

**Memory Usage:**
- 10,000 cells: 4-8 GB
- 100,000 cells: 16-32 GB
- 1,000,000 cells: 64+ GB

## Troubleshooting

### Issue: "Out of memory"

```bash
# Reduce number of cells or genes
# Or increase available RAM
# Or use sparse matrix formats
```

### Issue: "Poor UMAP separation"

```python
# Try different parameters
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50)
sc.tl.umap(adata, min_dist=0.1, spread=1.0)

# Or try different clustering resolution
sc.tl.leiden(adata, resolution=1.0)  # Higher = more clusters
```

### Issue: "Batch effects visible in UMAP"

```python
# Apply batch correction
import bbknn
bbknn.bbknn(adata, batch_key='sample')

# Or use Harmony
from harmony import harmonize
adata.obsm['X_harmony'] = harmonize(adata.X, adata.obs, batch_key='sample')
```

## References

- **Seurat:** Stuart et al. (2019) Cell
- **Scanpy:** Wolf et al. (2018) Genome Biology
- **Leiden Clustering:** Traag et al. (2019) Scientific Reports
- **UMAP:** McInnes et al. (2018) arXiv
- **Monocle3:** Cao et al. (2019) Nature Methods

## License

MIT License — See LICENSE file

## Contact

Adam Hoffman  
Email: adamhoffman21@hotmail.ca  
GitHub: [@adamhoffman2155-hue](https://github.com/adamhoffman2155-hue)

## Acknowledgments

- Test data from [10x Genomics](https://www.10xgenomics.com/)
- Marker gene databases: [PanglaoDB](https://panglaodb.se/), [ImmGen](http://www.immgen.org/)
- Workflow inspired by [Seurat tutorials](https://satijalab.org/seurat/) and [Scanpy documentation](https://scanpy.readthedocs.io/)

