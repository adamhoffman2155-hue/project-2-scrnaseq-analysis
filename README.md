# Project 2: scRNA-seq Tumour Microenvironment Analysis

[![CI](https://github.com/adamhoffman2155-hue/project-2-scrnaseq-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/adamhoffman2155-hue/project-2-scrnaseq-analysis/actions/workflows/ci.yml)

**Research question:** How does MSI status remodel the immune microenvironment in gastroesophageal adenocarcinoma?

This is the second project in a [computational biology portfolio](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio). After Project 1 surfaced immune-related pathways in bulk RNA-seq, this project digs into the tumour microenvironment at single-cell resolution to understand how MSI-high status reshapes immune infiltration — a question directly relevant to checkpoint inhibitor response in GEA.

## What It Does

`scripts/run_workflow.py` is a Python/Scanpy end-to-end orchestrator that runs nine steps on either a user-provided `.h5ad` or synthetic data:

1. **Load** — read `.h5ad` (or call `generate_synthetic_data()` when no input is passed)
2. **QC** — `qc_metrics.calculate_qc_metrics` + `filter_cells` (min 200 / max 5000 genes, ≤ 20 % MT, ≤ 50 % ribo) + `filter_genes` (≥ 3 cells)
3. **Normalize** — `sc.pp.normalize_total(target_sum=1e4)` + `log1p`
4. **HVG + PCA** — 2 000 HVGs (`seurat_v3`) + 50-component PCA
5. **Neighbors / UMAP / Leiden** — resolution configurable via `--resolution`
6. **Annotation** — marker-based (`annotation.annotate_by_markers`) over a default T/B/Monocyte/NK/DC panel + cluster-majority refinement
7. **Differential expression** — `sc.tl.rank_genes_groups` (Wilcoxon) by cluster
8. **Plots** — Leiden UMAP, cell-type UMAP, DE gene ranks (written under `<output>/plots/`)
9. **Save** — `processed.h5ad`, `cell_metadata.csv`, `de_results.csv`

Two standalone R scripts (`scripts/clustering.R` Seurat clustering, `scripts/trajectory.R` pseudotime) sit alongside the Python workflow for the specific steps where R tooling is preferred; they are not orchestrated by `run_workflow.py`.

## Methods & Tools

| Category | Tools |
|----------|-------|
| End-to-end workflow | Python + Scanpy (Leiden, UMAP, Wilcoxon DE) |
| Complementary R | Seurat clustering, pseudotime trajectory |
| Clustering | Graph-based / Leiden |
| Dimensionality reduction | PCA, UMAP |
| Visualisation | matplotlib, Scanpy plotting |
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
│   ├── run_workflow.py           # Python/Scanpy end-to-end orchestrator
│   ├── generate_synthetic_data.py# Reproducible synthetic input data
│   ├── qc_metrics.py             # QC metrics + filtering
│   ├── annotation.py             # Marker-based cell-type annotation
│   ├── visualization_utils.py    # Plot helpers
│   ├── clustering.R              # Standalone Seurat clustering
│   └── trajectory.R              # Standalone pseudotime trajectory
├── tests/                        # pytest suite
└── .gitignore
```

Run outputs (processed `.h5ad`, QC reports, figures, DE tables) are written under `data/` and `results/` at runtime and are gitignored.

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

# Synthetic end-to-end run (no input data required)
python scripts/run_workflow.py --output results/
# or with real data:
python scripts/run_workflow.py --input data/your.h5ad --output results/
```

## My Role

I scoped the MSI/TME question based on my thesis cohort, evaluated cluster assignments against known immune biology, and interpreted findings in the clinical context of GEA immunotherapy response. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 2 of 7**. It extends Project 1's bulk transcriptomics into single-cell resolution, asking whether the immune signatures detected in bulk are driven by specific cell populations. The immune subtype classifications from this work feed into the survival model in Project 6. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
