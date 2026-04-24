#!/usr/bin/env python3
"""
Main Executable Pipeline for scRNA-seq Analysis
Runs the full workflow: load -> QC -> normalize -> cluster -> annotate -> DE -> plot -> save
"""

import argparse
import os
import sys
from pathlib import Path

import scanpy as sc

# Ensure sibling modules are importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotation import annotate_by_markers, refine_annotations
from generate_synthetic_data import generate_synthetic_data
from qc_metrics import calculate_qc_metrics, filter_cells, filter_genes

# ---------------------------------------------------------------------------
# Default marker dictionary
# ---------------------------------------------------------------------------
DEFAULT_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "IL7R"],
    "B_cell": ["CD19", "MS4A1", "CD79A"],
    "Monocyte": ["CD14", "LYZ", "CST3"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1"],
    "Dendritic": ["FCER1A", "CLEC10A", "CD1C"],
}


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


def run_qc(adata):
    """Step 2: QC filtering."""
    print("[Step 2] Calculating QC metrics and filtering")
    adata = calculate_qc_metrics(adata)
    adata = filter_cells(adata, min_genes=200, max_genes=5000, max_mt=20, max_rb=50)
    adata = filter_genes(adata, min_cells=3)
    return adata


def normalize(adata):
    """Step 3: Normalize counts."""
    print("[Step 3] Normalizing (total-count + log1p)")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep raw for DE
    return adata


def select_hvg_and_pca(adata, n_top_genes=2000, n_pcs=50):
    """Step 4: HVG selection + PCA."""
    print(f"[Step 4] Selecting {n_top_genes} HVGs and running PCA ({n_pcs} PCs)")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer=None, subset=False
    )
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    return adata


def umap_and_cluster(adata, resolution=0.5):
    """Step 5: Neighbor graph, UMAP, and Leiden clustering."""
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


def differential_expression(adata, groupby="leiden"):
    """Step 7: DE analysis."""
    print(f"[Step 7] Differential expression (Wilcoxon, groupby={groupby})")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon")
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
    parser.add_argument(
        "--input", "-i", default=None, help="Path to input h5ad file (omit to use synthetic data)"
    )
    parser.add_argument(
        "--output", "-o", default="results", help="Output directory (default: results)"
    )
    parser.add_argument(
        "--config", "-c", default=None, help="Path to YAML config (optional, not yet implemented)"
    )
    parser.add_argument(
        "--resolution", type=float, default=0.5, help="Leiden clustering resolution (default: 0.5)"
    )
    parser.add_argument(
        "--n-hvgs", type=int, default=2000, help="Number of highly variable genes (default: 2000)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    adata = load_data(args.input)
    adata = run_qc(adata)
    adata = normalize(adata)
    adata = select_hvg_and_pca(adata, n_top_genes=args.n_hvgs)
    adata = umap_and_cluster(adata, resolution=args.resolution)
    adata = annotate(adata)
    adata = differential_expression(adata)
    generate_plots(adata, args.output)
    save_results(adata, args.output)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
