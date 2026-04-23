# Project 2: scRNA-seq Tumour Microenvironment Analysis

**Research question:** How does MSI status remodel the immune microenvironment in gastroesophageal adenocarcinoma?

This is the second project in a [computational biology portfolio](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio). After Project 1 surfaced immune-related pathways in bulk RNA-seq, this project digs into the tumour microenvironment at single-cell resolution to understand how MSI-high status reshapes immune infiltration — a question directly relevant to checkpoint inhibitor response in GEA.

## What It Does

Single-cell RNA-seq analysis of GEA tumour microenvironment:

1. **Quality control** — Cell filtering (nUMI, nGenes, %MT), doublet detection
2. **Normalization** — Library size normalization, log transformation, HVG selection
3. **Dimensionality reduction** — PCA, UMAP for visualization
4. **Clustering** — Leiden / graph-based clustering for cell population identification
5. **Cell type annotation** — Marker-based and reference-based annotation
6. **Differential expression** — Cluster-specific DE analysis
7. **Trajectory analysis** — Pseudotime ordering
8. **Visualisation utilities** — QC, UMAP, and marker-plot helpers

The workflow is orchestrated from `scripts/run_workflow.py` and is runnable on synthetic fixtures produced by `scripts/generate_synthetic_data.py`, so the pipeline can be exercised end-to-end without access to a private cohort.

## Methods & Tools

| Category | Tools |
|----------|-------|
| R analysis | Seurat (clustering, trajectory) |
| Python analysis | Scanpy-style QC + annotation utilities |
| Clustering | Graph-based / Leiden |
| Dimensionality reduction | PCA, UMAP |
| Visualisation | matplotlib, seaborn, ggplot2 |
| Testing | pytest |
| Environment | Docker, Conda |

## Project Structure

```
project-2-scrnaseq-analysis/
├── Dockerfile
├── environment.yml
├── config/
│   └── analysis_config.yaml
├── scripts/
│   ├── run_workflow.py           # End-to-end orchestrator
│   ├── generate_synthetic_data.py# Reproducible synthetic input data
│   ├── qc_metrics.py             # QC metrics (nUMI, nGenes, %MT, doublets)
│   ├── clustering.R              # Seurat clustering
│   ├── annotation.py             # Marker- and reference-based cell-type annotation
│   ├── trajectory.R              # Pseudotime trajectory
│   └── visualization_utils.py    # Plot helpers
├── tests/                        # pytest suite
└── .gitignore
```

Run outputs (QC reports, figures, tables, Seurat/AnnData objects) are written under `data/` and `results/` at runtime and are gitignored.

## Quick Start

```bash
git clone https://github.com/adamhoffman2155-hue/project-2-scrnaseq-analysis.git
cd project-2-scrnaseq-analysis

# Using Docker
docker build -t scrnaseq-analysis .
docker run -it -v $(pwd):/workspace scrnaseq-analysis bash

# Or Conda
conda env create -f environment.yml
conda activate scrnaseq-analysis

# Synthetic end-to-end run
python scripts/generate_synthetic_data.py
python scripts/run_workflow.py
```

## My Role

I scoped the MSI/TME question based on my thesis cohort, evaluated cluster assignments against known immune biology, and interpreted findings in the clinical context of GEA immunotherapy response. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 2 of 7**. It extends Project 1's bulk transcriptomics into single-cell resolution, asking whether the immune signatures detected in bulk are driven by specific cell populations. The immune subtype classifications from this work feed into the survival model in Project 6. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
