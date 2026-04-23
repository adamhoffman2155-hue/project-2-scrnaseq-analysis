# Project 2: scRNA-seq Tumour Microenvironment Analysis

> **A pipeline that looks at individual immune cells inside a tumor — one cell at a time — to understand why some tumors hide from the immune system.**

## The short version

**What this project does.** Takes single-cell RNA-sequencing data (where you get a gene-activity readout for every individual cell, not a bulk average) and sorts the cells into populations — T cells, B cells, monocytes, exhausted T cells, and so on — so you can ask how the cell makeup changes between tumor types.

**The question behind it.** Microsatellite-unstable (MSI-high) tumors accumulate more mutations, which in theory makes them more visible to the immune system — and that's why they sometimes respond to immunotherapy when other tumors don't. But "the tumor" is actually a mix of many cell types. To understand what MSI is really doing to the immune environment, you need to separate those cells.

**What the proof-of-concept shows.** Running the workflow on a standard 2,700-cell blood sample (the canonical benchmark for every single-cell tool), it correctly recovers the expected cell types — T cells, B cells, monocytes, NK cells — from gene-expression patterns alone. This is an **infrastructure check**, not the MSI/tumor biology answer.

**Why it matters.** Immunotherapy decisions increasingly depend on understanding the *composition* of the tumor microenvironment — how many T cells, which subtypes, whether they look exhausted — not just the cancer cells themselves. Single-cell tools make that possible at a resolution bulk RNA-seq can't reach.

---

_The rest of this README is technical detail for bioinformaticians, recruiters doing a deep review, or anyone reproducing the work._

## At a Glance

| | |
|---|---|
| **Stack** | Scanpy · Seurat v5 · Leiden · UMAP · Monocle3 · CellChat · Docker |
| **Data** | GEA tumour microenvironment (target); 10x PBMC 3k (POC — infrastructure check) |
| **POC headline** | 2,643 cells post-QC; 7 Leiden clusters; T/B/Monocyte/NK structure recovered |
| **Status** | Full analysis: **Full-data target** (needs tumor atlas); POC: **Runnable POC** |
| **Role** | MSI/TME framing; cluster-vs-biology review; implementation AI-assisted |
| **Portfolio** | Project 2 of 7 · [full narrative](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) |

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

## Methods & Tools

| Category | Tools |
|----------|-------|
| R Analysis | Seurat v5, Monocle3, CellChat |
| Python Analysis | Scanpy, scvi-tools |
| Clustering | Leiden algorithm |
| Visualization | UMAP, ggplot2, matplotlib, seaborn |
| Batch Correction | Harmony, BBKNN |
| Environment | Docker, Conda |

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

## Proof of Concept

A minimal end-to-end run of the scRNA-seq workflow on a small public dataset, to prove that the pipeline (QC → normalize → HVG → PCA → neighbors → UMAP → Leiden → marker annotation) actually runs on real 10x data and produces sensible biology.

**Scope — read this first:** this POC validates the workflow on standard PBMC data, **NOT** the MSI-GEA question. The full MSI-immune analysis requires a larger dataset (e.g. Pelka 2021 CRC atlas, ~2 GB) and is out of scope for a 1-day POC. Treat the POC as an "infrastructure check", not a biological result.

**Dataset:** 10x Genomics PBMC 3k (2,700 peripheral blood mononuclear cells from a healthy donor) — the canonical scanpy / Seurat tutorial dataset.

**Reproduce:**
```bash
pip install scanpy igraph leidenalg  # plus pandas numpy scipy matplotlib
python scripts/poc/run_poc.py
```

Outputs land in `results/poc/`:
- `umap_leiden.png` — UMAP coloured by Leiden cluster
- `umap_celltype.png` — UMAP coloured by marker-based cell type
- `cluster_composition.csv` — Leiden × cell-type count table
- `poc_summary.txt` — full run log

**Actual numbers from this run:**

| Metric | Value |
|---|---|
| N cells loaded | 2,700 |
| N cells after QC | 2,643 |
| N genes after QC | 13,714 |
| N Leiden clusters (res=0.5) | 7 |
| Runtime | ~18 s |

Marker-based cell-type distribution: T_cell (1,177), Monocyte (688), NK_cell (437), B_cell (341). Totals and major-lineage counts match the canonical PBMC3k tutorial within ~10–15%.

**Limitations:**
- PBMC3k is a healthy-donor benchmark, **not** a tumor sample.
- 5-type marker annotation does not distinguish CD4 vs CD8 T, classical vs non-classical monocytes, megakaryocytes, or cycling cells.
- The POC does not cover integration, batch correction, doublet detection, pseudotime, or cell-cell communication — those belong to the full analysis.

## My Role

I scoped the MSI/TME question based on my thesis cohort, evaluated cluster assignments against known immune biology, and interpreted findings in the clinical context of GEA immunotherapy response. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 2 of 7**. It extends Project 1's bulk transcriptomics into single-cell resolution, asking whether the immune signatures detected in bulk are driven by specific cell populations. The immune subtype classifications from this work feed into the survival model in Project 6. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
