#!/usr/bin/env python3
"""
Main Executable Pipeline for scRNA-seq Analysis
Runs the full workflow: load -> QC -> normalize -> cluster -> annotate -> DE -> plot -> save
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import yaml

# Ensure sibling modules are importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qc_metrics import calculate_qc_metrics, filter_cells, filter_genes
from annotation import annotate_by_markers, refine_annotations
from generate_synthetic_data import generate_synthetic_data


# ---------------------------------------------------------------------------
# Default marker dictionary
# ---------------------------------------------------------------------------
DEFAULT_MARKERS = {
    "T_cell":     ["CD3D", "CD3E", "IL7R"],
    "B_cell":     ["CD19", "MS4A1", "CD79A"],
    "Monocyte":   ["CD14", "LYZ", "CST3"],
    "NK_cell":    ["NKG7", "GNLY", "KLRD1"],
    "Dendritic":  ["FCER1A", "CLEC10A", "CD1C"],
}


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(path):
    """Load a YAML config file. Returns {} when path is None."""
    if not path:
        return {}
    cfg_path = Path(path)
    if not cfg_path.is_file():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")
    with cfg_path.open() as fh:
        return yaml.safe_load(fh) or {}


def _get(cfg, section, key, default):
    """Fetch a nested config value with a default fallback."""
    return cfg.get(section, {}).get(key, default) if cfg else default


def load_marker_file(path):
    """Load a marker-gene CSV with columns cell_type,gene. Returns a dict."""
    import csv
    markers = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ct = row.get("cell_type") or row.get("celltype")
            gene = row.get("gene")
            if not ct or not gene:
                continue
            markers.setdefault(ct, []).append(gene)
    return markers


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


def run_qc(adata, cfg=None):
    """Step 2: QC filtering."""
    print("[Step 2] Calculating QC metrics and filtering")
    adata = calculate_qc_metrics(adata)
    adata = filter_cells(
        adata,
        min_genes=_get(cfg, "qc", "min_genes", 200),
        max_genes=_get(cfg, "qc", "max_genes", 5000),
        max_mt=_get(cfg, "qc", "max_mt_pct", 20),
        max_rb=_get(cfg, "qc", "max_rb_pct", 50),
    )
    adata = filter_genes(
        adata,
        min_cells=_get(cfg, "qc", "min_cells_per_gene", 3),
    )
    return adata


def normalize(adata, cfg=None):
    """Step 3: Normalize counts."""
    target_sum = _get(cfg, "normalization", "target_sum", 1e4)
    log_transform = _get(cfg, "normalization", "log_transform", True)
    print(f"[Step 3] Normalizing (total-count target_sum={target_sum}"
          f"{', log1p' if log_transform else ''})")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    if log_transform:
        sc.pp.log1p(adata)
    adata.raw = adata  # keep raw for DE
    return adata


def select_hvg_and_pca(adata, n_top_genes=2000, n_pcs=50, cfg=None):
    """Step 4: HVG selection + PCA."""
    n_top_genes = _get(cfg, "feature_selection", "n_top_genes", n_top_genes)
    n_pcs = _get(cfg, "pca", "n_pcs", n_pcs)
    max_value = _get(cfg, "scaling", "max_value", 10)
    use_hv = _get(cfg, "pca", "use_highly_variable", True)
    print(f"[Step 4] Selecting {n_top_genes} HVGs and running PCA ({n_pcs} PCs)")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3",
                                layer=None, subset=False)
    sc.pp.scale(adata, max_value=max_value)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=use_hv)
    return adata


def umap_and_cluster(adata, resolution=0.5, cfg=None):
    """Step 5: Neighbor graph, UMAP, and Leiden clustering."""
    resolution = _get(cfg, "clustering", "resolution", resolution)
    n_neighbors = _get(cfg, "umap", "n_neighbors", 15)
    random_state = _get(cfg, "clustering", "random_state", 42)
    method = _get(cfg, "clustering", "method", "leiden")
    print(f"[Step 5] Building neighbor graph, UMAP, and {method} (res={resolution})")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=30)
    sc.tl.umap(adata)
    if method == "louvain":
        sc.tl.louvain(adata, resolution=resolution, random_state=random_state)
        adata.obs["leiden"] = adata.obs["louvain"]  # downstream expects 'leiden'
    else:
        sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
    print(f"  Found {adata.obs['leiden'].nunique()} clusters")
    return adata


def annotate(adata, marker_dict=None, cfg=None):
    """Step 6: Cell type annotation."""
    print("[Step 6] Annotating cell types by marker genes")
    markers = marker_dict
    if markers is None:
        marker_file = _get(cfg, "annotation", "marker_file", None)
        if marker_file and Path(marker_file).is_file():
            print(f"  Loading markers from {marker_file}")
            markers = load_marker_file(marker_file)
        else:
            if marker_file:
                print(f"  Marker file {marker_file} not found; using DEFAULT_MARKERS")
            markers = DEFAULT_MARKERS
    adata = annotate_by_markers(adata, markers)
    adata = refine_annotations(adata, cluster_key="leiden", label_key="cell_type")
    ct_counts = adata.obs["cell_type_refined"].value_counts()
    for ct, n in ct_counts.items():
        print(f"  {ct}: {n} cells")
    return adata


def differential_expression(adata, groupby="leiden", cfg=None):
    """Step 7: DE analysis."""
    method = _get(cfg, "differential_expression", "method", "wilcoxon")
    print(f"[Step 7] Differential expression ({method}, groupby={groupby})")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)
    return adata


def generate_plots(adata, output_dir):
    """Step 8: Generate all plots."""
    print("[Step 8] Generating plots")
    plot_dir = Path(output_dir) / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(plot_dir)

    sc.pl.umap(adata, color=["leiden"], show=False,
               save="_leiden.png")
    if "cell_type_refined" in adata.obs.columns:
        sc.pl.umap(adata, color=["cell_type_refined"], show=False,
                   save="_celltype.png")
    sc.pl.rank_genes_groups(adata, n_genes=10, show=False,
                            save="_de_genes.png")
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
    parser = argparse.ArgumentParser(
        description="scRNA-seq analysis pipeline"
    )
    parser.add_argument("--input", "-i", default=None,
                        help="Path to input h5ad file (omit to use synthetic data)")
    parser.add_argument("--output", "-o", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--config", "-c", default=None,
                        help="Path to YAML config (e.g. config/analysis_config.yaml). "
                             "When provided, overrides --resolution and --n-hvgs defaults.")
    parser.add_argument("--resolution", type=float, default=0.5,
                        help="Leiden clustering resolution (default: 0.5)")
    parser.add_argument("--n-hvgs", type=int, default=2000,
                        help="Number of highly variable genes (default: 2000)")
    return parser.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    if cfg:
        print(f"[Config] Loaded {args.config}")

    adata = load_data(args.input)
    adata = run_qc(adata, cfg=cfg)
    adata = normalize(adata, cfg=cfg)
    adata = select_hvg_and_pca(adata, n_top_genes=args.n_hvgs, cfg=cfg)
    adata = umap_and_cluster(adata, resolution=args.resolution, cfg=cfg)
    adata = annotate(adata, cfg=cfg)
    adata = differential_expression(adata, cfg=cfg)
    generate_plots(adata, args.output)
    save_results(adata, args.output)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
