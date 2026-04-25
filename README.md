# Project 2: scRNA-seq Tumour Microenvironment Analysis

![CI](https://github.com/adamhoffman2155-hue/project-2-scrnaseq-analysis/actions/workflows/ci.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.11-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Repro](https://img.shields.io/badge/FAIR_DOME_CURE-11%2F14_%7C_5%2F7_%7C_4%2F4-brightgreen)

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
├── README.md
├── Dockerfile
├── environment.yml
├── config/
│   └── analysis_config.yaml
├── scripts/
│   ├── run_workflow.py                  # Python/Scanpy end-to-end orchestrator
│   ├── generate_synthetic_data.py       # Reproducible synthetic input data
│   ├── qc_metrics.py                    # QC metrics + filtering
│   ├── annotation.py                    # Marker-based cell-type annotation
│   ├── visualization_utils.py           # Plot helpers
│   ├── benchmark_harmony_vs_scvi.py     # Harmony vs scVI batch-correction benchmark
│   ├── clustering.R                     # Standalone Seurat clustering
│   ├── trajectory.R                     # Standalone pseudotime trajectory
│   └── poc/
│       └── run_poc.py                   # Proof-of-concept runner
├── tests/                               # pytest suite
└── results/
    └── poc/                             # POC outputs (committed)
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

## Proof of Concept

A minimal end-to-end run of the scRNA-seq workflow on a small public dataset, to prove that the pipeline (QC → normalize → HVG → PCA → neighbors → UMAP → Leiden → marker annotation) actually runs on real 10x data and produces sensible biology.

**Scope — read this first:** this POC validates the workflow on standard PBMC data, **NOT** the MSI-GEA question. The full MSI-immune analysis requires a larger dataset (e.g. Pelka 2021 CRC atlas, ~2 GB) and is out of scope for a 1-day POC. Treat the POC as an "infrastructure check", not a biological result.

**Dataset:** 10x Genomics PBMC 3k (2,700 peripheral blood mononuclear cells from a healthy donor) — the canonical scanpy / Seurat tutorial dataset.
Upstream URL: <https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>

The script tries `sc.datasets.pbmc3k()` first and falls back to a public GitHub mirror of the raw 10x `barcodes.tsv` / `genes.tsv` / `matrix.mtx` files if the canonical URL is unreachable from the execution environment.

**Reproduce:**

```bash
pip install scanpy igraph leidenalg  # plus pandas numpy scipy matplotlib
python scripts/poc/run_poc.py
```

Outputs land in `results/poc/`:

- `umap_leiden.png` — UMAP coloured by Leiden cluster
- `umap_celltype.png` — UMAP coloured by marker-based cell type
- `cluster_composition.csv` — Leiden × cell-type count table
- `poc_summary.txt` — full run log (N cells, N clusters, top markers, caveats)

**Actual numbers from this run:**

| Metric | Value |
|---|---|
| N cells loaded | 2,700 |
| N genes loaded | 32,738 |
| N cells after QC (`min_genes=200`, `pct_mt<5`) | 2,643 |
| N genes after QC (`min_cells=3`) | 13,714 |
| N Leiden clusters (res=0.5) | 7 |
| Runtime | ~18 s |

Marker-based cell-type distribution:

| Cell type | N cells |
|---|---|
| T_cell | 1,177 |
| Monocyte | 688 |
| NK_cell | 437 |
| B_cell | 341 |

**Cluster → cell type (top-5 t-test markers):**

| Cluster | N | Assigned | Top markers |
|---|---|---|---|
| 0 | 1177 | T_cell | LDHB, CD3D, RPS12, RPS27, RPS25 |
| 1 |  341 | B_cell | CD74, HLA-DRA, CD79A, HLA-DPB1, HLA-DRB1 |
| 2 |  640 | Monocyte | FTL, CST3, TYROBP, FTH1, AIF1 |
| 3 |   12 | Monocyte* | SDPR, GPX1, TAGLN2, GNG11, PF4 |
| 4 |  429 | NK_cell | NKG7, B2M, GZMA, CST7, CCL5 |
| 5 |   36 | Monocyte | HLA-DRA, CD74, CST3, HLA-DPA1, HLA-DPB1 |
| 6 |    8 | NK_cell* | DUT, PSME2, GAPDH, CFL1, ITGB1BP1 |

*Labels for small clusters are known to be wrong — see "Limitations".*

Totals and major-lineage counts match the canonical PBMC3k tutorial within ~10–15%. Minor subpopulations (dendritic cells, megakaryocytes) are present as distinct small clusters in the data but not distinguishable with this 5-type marker panel.

**Limitations:**

- PBMC3k is a healthy-donor benchmark, **not** a tumor sample. There is no MSI/tumor metadata in this dataset, and no immune-exhaustion / checkpoint biology to recover.
- Cell-type annotation uses a 5-type argmax over mean marker expression. It does not distinguish CD4 vs CD8 T, classical vs non-classical monocytes, megakaryocytes (no PF4/PPBP marker set), or cycling cells. Small clusters (cluster 3: megakaryocyte, cluster 6: likely cycling) get mis-labelled by the argmax rule — inspect `rank_genes_groups` before trusting the label column.
- Leiden at `resolution=0.5` returns 7 clusters vs the tutorial's 8; raising resolution would split CD4/CD8 T.
- The POC does not cover integration, batch correction, doublet detection, pseudotime, or cell-cell communication — those belong to the full analysis.

## My Role

I scoped the MSI/TME question based on my thesis cohort, evaluated cluster assignments against known immune biology, and interpreted findings in the clinical context of GEA immunotherapy response. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 2 of 7**. It extends Project 1's bulk transcriptomics into single-cell resolution, asking whether the immune signatures detected in bulk are driven by specific cell populations. The immune subtype classifications from this work feed into the survival model in Project 6. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

### Cross-project data flow

```
Project 1 (bulk RNA-seq DE)     ──┐  immune signatures
                                  │
Project 2 (this one — scRNA-seq) ─┼──▶  Project 6 (survival — TME subtype covariate)
                                  │
```

- **Upstream** — consumes sparse count matrices (`.h5ad`) from TCGA-STAD or user input; synthetic data generator available offline.
- **Downstream** — TME immune-subtype labels feed the Cox covariate panel in Project 6 (narrative input).

## Benchmarks

| Benchmark | Output | Summary |
| --- | --- | --- |
| Harmony vs scVI batch correction | [`results/benchmark/batch_correction_metrics.md`](results/benchmark/batch_correction_metrics.md) | On a 600-cell synthetic AnnData with 2 simulated batches, Harmony compresses batch silhouette toward 0 while preserving cell-type structure. scVI is guarded (optional `scvi-tools` install) so CI passes without torch. |

Rebuild with `python scripts/benchmark_harmony_vs_scvi.py`.

## Reproducibility

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the FAIR-BioRS / DOME / CURE self-scorecard (11/14 · 5/7 · 4/4).

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
