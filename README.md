# Project 2: scRNA-seq Tumour Microenvironment Analysis

**Research question:** How does MSI status remodel the immune microenvironment in gastroesophageal adenocarcinoma?

This is the second project in a [computational biology portfolio](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio). After Project 1 surfaced immune-related pathways in bulk RNA-seq, this project digs into the tumour microenvironment at single-cell resolution to understand how MSI-high status reshapes immune infiltration — a question directly relevant to checkpoint inhibitor response in GEA.

## What It Does

Single-cell RNA-seq analysis of GEA tumour microenvironment:

1. **Quality control** — Cell filtering (nUMI, nGenes, %MT), doublet detection
2. **Normalization** — Library size normalization, log transformation, HVG selection
3. **Dimensionality reduction** — PCA, UMAP for visualization
4. **Clustering** — Leiden algorithm for cell population identification
5. **Cell type annotation** — Marker-based and reference-based annotation
6. **Differential expression** — Cluster-specific DE analysis
7. **Trajectory analysis** — Pseudotime ordering with Monocle3
8. **Cell-cell communication** — CellChat interaction analysis

The analysis explores MSI-high immune infiltration patterns and T-cell exhaustion signatures relevant to checkpoint inhibitor response in GEA.

## Methods & Tools

| Category | Tools |
|----------|-------|
| R Analysis | Seurat v5, Monocle3, CellChat |
| Python Analysis | Scanpy, scvi-tools |
| Clustering | Leiden algorithm |
| Visualization | UMAP, ggplot2, matplotlib, seaborn |
| Batch Correction | Harmony, BBKNN |
| Environment | Docker, Conda |

## Project Structure

```
project-2-scrnaseq-analysis/
├── README.md
├── Dockerfile
├── environment.yml
├── config/
│   └── analysis_config.yaml
├── scripts/
│   ├── annotation.py
│   ├── clustering.R
│   ├── generate_synthetic_data.py
│   ├── qc_metrics.py
│   ├── run_workflow.py
│   ├── trajectory.R
│   └── visualization_utils.py
└── tests/
    ├── __init__.py
    └── test_workflow.py
```

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

# Run the complete workflow
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
