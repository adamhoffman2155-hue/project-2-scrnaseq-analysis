"""
Proof-of-Concept: scRNA-seq workflow on 10x Genomics PBMC 3k.

This POC validates the scRNA-seq pipeline (QC -> normalize -> HVG -> PCA ->
neighbors -> UMAP -> Leiden -> marker annotation) end-to-end on real 10x data.

IMPORTANT: This POC does NOT answer the MSI / tumor-immune microenvironment
biological question. It is an "infrastructure POC" that proves the pipeline
runs on real 10x data and produces canonical PBMC cell-type structure. The
full MSI-GEA analysis requires a much larger dataset (e.g. Pelka 2021 CRC
atlas, ~2 GB) which is out of scope for a 1-day POC.

Run:
    python scripts/poc/run_poc.py
"""

from __future__ import annotations

import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results" / "poc"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Local 10x mtx directory (hg19). Populated via download_pbmc3k() if missing.
DATA_DIR = REPO_ROOT / "data" / "raw" / "hg19"

SUMMARY_PATH = RESULTS_DIR / "poc_summary.txt"
UMAP_LEIDEN_PATH = RESULTS_DIR / "umap_leiden.png"
UMAP_CELLTYPE_PATH = RESULTS_DIR / "umap_celltype.png"
COMPOSITION_PATH = RESULTS_DIR / "cluster_composition.csv"


# ---------------------------------------------------------------------------
# Scanpy config
# ---------------------------------------------------------------------------
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor="white")

# Canonical PBMC marker genes
MARKERS = {
    "T_cell":    ["CD3D", "CD3E", "IL7R"],
    "B_cell":    ["CD19", "MS4A1", "CD79A"],
    "Monocyte":  ["CD14", "LYZ", "CST3"],
    "NK_cell":   ["NKG7", "GNLY", "KLRD1"],
    "Dendritic": ["FCER1A", "CLEC10A"],
}


def log(msg: str) -> None:
    print(f"[poc] {msg}", flush=True)


# Fallback raw 10x mtx mirror (read-only, public GitHub copy of the
# canonical PBMC3k filtered_gene_bc_matrices/hg19 files from 10x Genomics).
# Original upstream: https://cf.10xgenomics.com/samples/cell-exp/1.1.0/
#   pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
PBMC3K_MIRROR = (
    "https://raw.githubusercontent.com/outobi/pbmc_Seurat_scRNAseq/"
    "c3966e791982c3b1755a361f28c8080752f2337d/hg19"
)


def download_pbmc3k_mtx(target_dir: Path) -> None:
    """Download PBMC3k 10x mtx files (barcodes/genes/matrix) to target_dir."""
    import urllib.request

    target_dir.mkdir(parents=True, exist_ok=True)
    for fname in ("barcodes.tsv", "genes.tsv", "matrix.mtx"):
        out = target_dir / fname
        if out.exists() and out.stat().st_size > 0:
            continue
        url = f"{PBMC3K_MIRROR}/{fname}"
        log(f"Downloading {fname} from {url}")
        urllib.request.urlretrieve(url, out)


def load_pbmc3k() -> "sc.AnnData":
    """
    Load PBMC3k. Tries sc.datasets.pbmc3k() first (which hits scanpy's
    canonical URL); if that is blocked by a sandbox, falls back to a
    public GitHub mirror of the raw 10x mtx files.
    """
    try:
        log("Trying sc.datasets.pbmc3k() ...")
        adata = sc.datasets.pbmc3k()
        return adata
    except Exception as e:
        log(f"sc.datasets.pbmc3k() failed: {e!r}")
        log("Falling back to raw 10x mtx mirror on GitHub.")
        download_pbmc3k_mtx(DATA_DIR)
        adata = sc.read_10x_mtx(
            str(DATA_DIR), var_names="gene_symbols", cache=False
        )
        return adata


def main() -> int:
    t0 = time.time()
    caveats: list[str] = []

    # -----------------------------------------------------------------
    # 1. Load PBMC3k
    # -----------------------------------------------------------------
    adata = load_pbmc3k()
    adata.var_names_make_unique()
    n_cells_loaded = adata.n_obs
    n_genes_loaded = adata.n_vars
    log(f"Loaded: {n_cells_loaded} cells x {n_genes_loaded} genes")

    # -----------------------------------------------------------------
    # 2. QC metrics
    # -----------------------------------------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # -----------------------------------------------------------------
    # 3. Filter cells and genes
    # -----------------------------------------------------------------
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < 5, :].copy()
    n_cells_post_qc = adata.n_obs
    n_genes_post_qc = adata.n_vars
    log(f"Post-QC: {n_cells_post_qc} cells x {n_genes_post_qc} genes")

    # -----------------------------------------------------------------
    # 4. Normalize, HVG, scale
    # -----------------------------------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Keep a raw copy for marker gene plotting
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    # -----------------------------------------------------------------
    # 5. PCA / neighbors / UMAP / Leiden
    # -----------------------------------------------------------------
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata, resolution=0.5, flavor="igraph",
                     n_iterations=2, directed=False)
    except Exception:
        # Fall back to default leiden (leidenalg) if igraph flavor unavailable
        sc.tl.leiden(adata, resolution=0.5)
    n_clusters = adata.obs["leiden"].nunique()
    log(f"Leiden clusters: {n_clusters}")

    # -----------------------------------------------------------------
    # 6. Marker-based cell-type annotation
    # -----------------------------------------------------------------
    # Score each cluster by mean expression of each marker set,
    # then assign the cluster to the highest-scoring type.
    raw = adata.raw.to_adata() if adata.raw is not None else adata

    cluster_scores: dict[str, dict[str, float]] = {}
    for cl in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        mask = (adata.obs["leiden"] == cl).values
        sub = raw[mask, :]
        scores: dict[str, float] = {}
        for ctype, genes in MARKERS.items():
            present = [g for g in genes if g in sub.var_names]
            if not present:
                scores[ctype] = float("nan")
                continue
            X = sub[:, present].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            scores[ctype] = float(np.mean(X))
        cluster_scores[cl] = scores

    scores_df = pd.DataFrame(cluster_scores).T  # rows = clusters
    scores_df.index.name = "leiden"
    cluster_to_celltype = {
        cl: scores_df.loc[cl].idxmax() for cl in scores_df.index
    }
    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(cluster_to_celltype).astype("category")
    )

    # Rank genes per cluster for the summary (top markers)
    sc.tl.rank_genes_groups(adata, "leiden", method="t-test", n_genes=5)
    top_markers = {}
    names = adata.uns["rank_genes_groups"]["names"]
    for cl in names.dtype.names:
        top_markers[cl] = list(names[cl])[:5]

    # -----------------------------------------------------------------
    # 7. Plots
    # -----------------------------------------------------------------
    sc.pl.umap(adata, color="leiden", legend_loc="on data",
               title="PBMC3k Leiden (res=0.5)", show=False)
    plt.savefig(UMAP_LEIDEN_PATH, bbox_inches="tight", dpi=150)
    plt.close("all")

    sc.pl.umap(adata, color="cell_type",
               title="PBMC3k marker-based cell types", show=False)
    plt.savefig(UMAP_CELLTYPE_PATH, bbox_inches="tight", dpi=150)
    plt.close("all")

    # -----------------------------------------------------------------
    # 8. Cluster composition (cluster x cell_type counts)
    # -----------------------------------------------------------------
    comp = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])
    comp.to_csv(COMPOSITION_PATH)

    celltype_counts = adata.obs["cell_type"].value_counts().to_dict()

    # -----------------------------------------------------------------
    # 9. Summary file
    # -----------------------------------------------------------------
    elapsed = time.time() - t0
    lines: list[str] = []
    lines.append("PROOF OF CONCEPT - scRNA-seq workflow on PBMC3k")
    lines.append("=" * 60)
    lines.append("")
    lines.append("Scope note:")
    lines.append("  This POC validates the scRNA-seq infrastructure only.")
    lines.append("  It does NOT answer the MSI / tumor-immune question.")
    lines.append("  The full MSI-GEA analysis needs a larger dataset")
    lines.append("  (e.g. Pelka 2021 CRC atlas, ~2 GB) and is out of")
    lines.append("  scope for a 1-day POC.")
    lines.append("")
    lines.append("Dataset: 10x Genomics PBMC 3k")
    lines.append(
        "  Upstream: https://cf.10xgenomics.com/samples/cell-exp/1.1.0/"
        "pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
    )
    lines.append(
        "  Loader: sc.datasets.pbmc3k() with fallback to raw 10x mtx"
        " files (barcodes/genes/matrix) mirrored on GitHub."
    )
    lines.append("")
    lines.append("Counts:")
    lines.append(f"  N cells loaded         : {n_cells_loaded}")
    lines.append(f"  N genes loaded         : {n_genes_loaded}")
    lines.append(f"  N cells after QC       : {n_cells_post_qc}")
    lines.append(f"  N genes after QC       : {n_genes_post_qc}")
    lines.append(f"  N Leiden clusters      : {n_clusters}")
    lines.append("")
    lines.append("Cell-type distribution (marker-based):")
    for ctype, n in sorted(celltype_counts.items(), key=lambda kv: -kv[1]):
        lines.append(f"  {ctype:<12s} : {n}")
    lines.append("")
    lines.append("Cluster -> cell type:")
    for cl, ctype in sorted(
        cluster_to_celltype.items(), key=lambda kv: int(kv[0])
    ):
        n = int((adata.obs["leiden"] == cl).sum())
        lines.append(f"  cluster {cl:>2s} (n={n:>4d}) -> {ctype}")
    lines.append("")
    lines.append("Top 5 markers per Leiden cluster (t-test):")
    for cl in sorted(top_markers.keys(), key=lambda x: int(x)):
        lines.append(f"  cluster {cl:>2s}: {', '.join(top_markers[cl])}")
    lines.append("")
    lines.append("Caveats:")
    lines.append("  - PBMC3k is a benchmarking dataset, not a tumor sample.")
    lines.append("  - No MSI / tumor metadata available here.")
    lines.append("  - Marker-based annotation uses mean expression of a")
    lines.append("    small canonical gene set; borderline clusters may")
    lines.append("    collapse onto the dominant type (e.g. CD4/CD8 T merge).")
    lines.append("  - The canonical PBMC3k has 8 clusters including a")
    lines.append("    tiny megakaryocyte (PF4/PPBP) cluster and a")
    lines.append("    dendritic (FCER1A/CLEC10A) cluster. With this marker")
    lines.append("    set (no MK, small DC set) those clusters are")
    lines.append("    mis-assigned by argmax — e.g. in one run cluster 3")
    lines.append("    (SDPR/PF4/GNG11/PPBP markers = megakaryocyte) was")
    lines.append("    labelled Monocyte. Inspect rank_genes_groups output")
    lines.append("    before trusting the label column for small clusters.")
    for c in caveats:
        lines.append(f"  - {c}")
    lines.append("")
    lines.append(f"Elapsed: {elapsed:.1f} s")

    SUMMARY_PATH.write_text("\n".join(lines) + "\n")
    log(f"Wrote summary: {SUMMARY_PATH}")
    log(f"Done in {elapsed:.1f} s")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
