#!/usr/bin/env python3
"""
Main Executable Pipeline for scRNA-seq Analysis

Runs the full workflow: load -> QC -> normalize -> cluster -> annotate -> DE -> plot -> save.
Parameters can be overridden by a YAML config file via --config.
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

# Ensure sibling modules are importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qc_metrics import calculate_qc_metrics, filter_cells, filter_genes
from annotation import annotate_by_markers, refine_annotations
from generate_synthetic_data import generate_synthetic_data


# ---------------------------------------------------------------------------
# Default marker dictionary (used when config doesn't provide one)
# ---------------------------------------------------------------------------
DEFAULT_MARKERS = {
    "T_cell":     ["CD3D", "CD3E", "IL7R"],
    "B_cell":     ["CD19", "MS4A1", "CD79A"],
    "Monocyte":   ["CD14", "LYZ", "CST3"],
    "NK_cell":    ["NKG7", "GNLY", "KLRD1"],
    "Dendritic":  ["FCER1A", "CLEC10A", "CD1C"],
}

DEFAULT_CONFIG = {
    "qc": {
        "min_genes": 200,
        "max_genes": 5000,
        "max_mt_pct": 20,
        "max_rb_pct": 50,
        "min_cells_per_gene": 3,
    },
    "feature_selection": {"n_top_genes": 2000},
    "pca": {"n_pcs": 50},
    "clustering": {"resolution": 0.5},
    "differential_expression": {"method": "wilcoxon"},
}


def load_config(config_path):
    """Load YAML config and merge with defaults. Returns nested dict."""
    if not config_path:
        return DEFAULT_CONFIG

    try:
        import yaml
    except ImportError:
        print("Warning: pyyaml not installed; using default config", file=sys.stderr)
        return DEFAULT_CONFIG

    if not os.path.isfile(config_path):
        print(f"Warning: config file not found ({config_path}); using defaults", file=sys.stderr)
        return DEFAULT_CONFIG

    with open(config_path, "r") as fh:
        user_cfg = yaml.safe_load(fh) or {}

    # Shallow-merge top-level sections onto defaults
    merged = {k: dict(v) for k, v in DEFAULT_CONFIG.items()}
    for section, values in user_cfg.items():
        if isinstance(values, dict):
            merged.setdefault(section, {}).update(values)
        else:
            merged[section] = values
    return merged


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def load_data(input_path):
    """Step 1: Load data from h5ad or generate synthetic data."""
    if input_path and os.path.isfile(input_path):
        print(f"[Step 1] Loading data from {input_path}")
        adata = sc.read_h5ad(input_path)
    else:
        print("[Step 1] No input file found. Generating synthetic dataset.")
        adata = generate_synthetic_data()
    print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def run_qc(adata, config):
    """Step 2: QC filtering using config thresholds."""
    qc_cfg = config.get("qc", {})
    print("[Step 2] Calculating QC metrics and filtering")
    adata = calculate_qc_metrics(adata)
    adata = filter_cells(
        adata,
        min_genes=qc_cfg.get("min_genes", 200),
        max_genes=qc_cfg.get("max_genes", 5000),
        max_mt=qc_cfg.get("max_mt_pct", 20),
        max_rb=qc_cfg.get("max_rb_pct", 50),
    )
    adata = filter_genes(adata, min_cells=qc_cfg.get("min_cells_per_gene", 3))
    return adata


def normalize(adata):
    """Step 3: Normalize counts."""
    print("[Step 3] Normalizing (total-count + log1p)")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep raw for DE
    return adata


def select_hvg_and_pca(adata, config):
    """Step 4: HVG selection + PCA."""
    n_top_genes = config.get("feature_selection", {}).get("n_top_genes", 2000)
    n_pcs = config.get("pca", {}).get("n_pcs", 50)
    print(f"[Step 4] Selecting {n_top_genes} HVGs and running PCA ({n_pcs} PCs)")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3",
                                layer=None, subset=False)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    return adata


def umap_and_cluster(adata, config):
    """Step 5: Neighbor graph, UMAP, and Leiden clustering."""
    resolution = config.get("clustering", {}).get("resolution", 0.5)
    print(f"[Step 5] Building neighbor graph, UMAP, and Leiden (res={resolution})")
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)
    print(f"  Found {adata.obs['leiden'].nunique()} clusters")
    return adata


def annotate(adata, marker_dict=None):
    """Step 6: Cell type annotation."""
    print("[Step 6] Annotating cell types by marker genes")
    markers = marker_dict or DEFAULT_MARKERS
    adata = annotate_by_markers(adata, markers)
    adata = refine_annotations(adata, cluster_key="leiden", label_key="cell_type")
    ct_counts = adata.obs["cell_type_refined"].value_counts()
    for ct, n in ct_counts.items():
        print(f"  {ct}: {n} cells")
    return adata


def differential_expression(adata, config):
    """Step 7: DE analysis."""
    method = config.get("differential_expression", {}).get("method", "wilcoxon")
    print(f"[Step 7] Differential expression ({method}, groupby=leiden)")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method=method)
    return adata


def generate_plots(adata, output_dir):
    """Step 8: Generate all plots."""
    print("[Step 8] Generating plots")
    plot_dir = Path(output_dir) / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(plot_dir)

    sc.pl.umap(adata, color=["leiden"], show=False, save="_leiden.png")
    if "cell_type_refined" in adata.obs.columns:
        sc.pl.umap(adata, color=["cell_type_refined"], show=False, save="_celltype.png")
    sc.pl.rank_genes_groups(adata, n_genes=10, show=False, save="_de_genes.png")
    print(f"  Plots saved to {plot_dir}")


def save_results(adata, output_dir):
    """Step 9: Save processed data and tables."""
    print("[Step 9] Saving results")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    adata.write_h5ad(out / "processed.h5ad")
    adata.obs.to_csv(out / "cell_metadata.csv")

    if "rank_genes_groups" in adata.uns:
        de_df = sc.get.rank_genes_groups_df(adata, group=None)
        de_df.to_csv(out / "de_results.csv", index=False)

    print(f"  Results saved to {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="scRNA-seq analysis pipeline")
    parser.add_argument("--input", "-i", default=None,
                        help="Path to input h5ad file (omit to use synthetic data)")
    parser.add_argument("--output", "-o", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--config", "-c", default="config/analysis_config.yaml",
                        help="Path to YAML config (default: config/analysis_config.yaml)")
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)

    adata = load_data(args.input)
    adata = run_qc(adata, config)
    adata = normalize(adata)
    adata = select_hvg_and_pca(adata, config)
    adata = umap_and_cluster(adata, config)
    adata = annotate(adata)
    adata = differential_expression(adata, config)
    generate_plots(adata, args.output)
    save_results(adata, args.output)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
