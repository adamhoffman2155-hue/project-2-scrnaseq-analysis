"""
Plotting Utilities for Single-Cell RNA-seq Analysis
Provides UMAP, violin, heatmap, and QC summary plots using scanpy and matplotlib.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _load_adata(adata_or_path):
    """Return an AnnData object, loading from path if necessary."""
    if isinstance(adata_or_path, (str, Path)):
        return sc.read_h5ad(str(adata_or_path))
    return adata_or_path


def _save_figure(fig, save_path):
    """Save a matplotlib figure if a path is provided."""
    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved figure to {save_path}")


# ---------------------------------------------------------------------------
# UMAP
# ---------------------------------------------------------------------------


def plot_umap(adata_or_path, color_by="leiden", save_path=None, **kwargs):
    """
    Plot UMAP colored by a given observation key or gene name.

    Parameters
    ----------
    adata_or_path : AnnData or str/Path
        Annotated data matrix or path to an h5ad file.
    color_by : str or list of str
        Key(s) in adata.obs or gene names to color by.
    save_path : str or None
        Path to save the figure.
    **kwargs
        Extra keyword arguments forwarded to ``sc.pl.umap``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    adata = _load_adata(adata_or_path)

    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP coordinates not found. Run sc.tl.umap first.")

    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))
    sc.pl.umap(adata, color=color_by, ax=ax, show=False, **kwargs)
    ax.set_title(f"UMAP colored by {color_by}")
    _save_figure(fig, save_path)
    return fig


# ---------------------------------------------------------------------------
# Violin
# ---------------------------------------------------------------------------


def plot_violin(adata_or_path, genes, groupby="leiden", save_path=None, **kwargs):
    """
    Violin plot of gene expression across groups.

    Parameters
    ----------
    adata_or_path : AnnData or str/Path
    genes : list of str
        Gene names to plot.
    groupby : str
        Observation key to group violins by.
    save_path : str or None
    **kwargs
        Extra keyword arguments forwarded to ``sc.pl.violin``.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    adata = _load_adata(adata_or_path)

    n_genes = len(genes)
    fig, axes = plt.subplots(1, n_genes, figsize=(4 * n_genes, 4))
    if n_genes == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes, strict=False):
        if gene not in adata.var_names:
            ax.set_title(f"{gene} (not found)")
            continue
        sc.pl.violin(adata, keys=gene, groupby=groupby, ax=ax, show=False, **kwargs)
        ax.set_title(gene)

    fig.tight_layout()
    _save_figure(fig, save_path)
    return fig


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------


def plot_heatmap(adata_or_path, markers, groupby="leiden", save_path=None, **kwargs):
    """
    Heatmap of marker gene expression per group (scanpy.pl style).

    Parameters
    ----------
    adata_or_path : AnnData or str/Path
    markers : list of str or dict
        Gene names. If dict, keys are group labels and values are gene lists.
    groupby : str
    save_path : str or None
    **kwargs
        Extra keyword arguments forwarded to ``sc.pl.heatmap``.

    Returns
    -------
    result : dict
        The dict returned by ``sc.pl.heatmap`` (contains axes handles).
    """
    adata = _load_adata(adata_or_path)

    if isinstance(markers, dict):
        var_names = markers
    else:
        present = [g for g in markers if g in adata.var_names]
        if not present:
            raise ValueError("None of the requested marker genes found in adata.var_names")
        var_names = present

    result = sc.pl.heatmap(
        adata,
        var_names=var_names,
        groupby=groupby,
        show=False,
        **kwargs,
    )

    if save_path is not None:
        fig = plt.gcf()
        _save_figure(fig, save_path)

    return result


# ---------------------------------------------------------------------------
# QC summary
# ---------------------------------------------------------------------------


def plot_qc_summary(adata_or_path, save_path=None):
    """
    QC summary panel: nUMI, nGenes, and pct_mt distributions.

    Parameters
    ----------
    adata_or_path : AnnData or str/Path
    save_path : str or None

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    adata = _load_adata(adata_or_path)

    metrics = []
    labels = []
    if "total_counts" in adata.obs.columns:
        metrics.append("total_counts")
        labels.append("nUMI (total_counts)")
    if "n_genes_by_counts" in adata.obs.columns:
        metrics.append("n_genes_by_counts")
        labels.append("nGenes")
    if "pct_mt" in adata.obs.columns:
        metrics.append("pct_mt")
        labels.append("% Mitochondrial")

    n = len(metrics)
    if n == 0:
        raise ValueError(
            "No QC metrics found in adata.obs. "
            "Run sc.pp.calculate_qc_metrics or qc_metrics.calculate_qc_metrics first."
        )

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]

    for ax, metric, label in zip(axes, metrics, labels, strict=False):
        values = adata.obs[metric].values
        ax.hist(values, bins=50, edgecolor="black", alpha=0.7)
        ax.axvline(
            np.median(values), color="red", linestyle="--", label=f"median={np.median(values):.1f}"
        )
        ax.set_xlabel(label)
        ax.set_ylabel("Frequency")
        ax.set_title(label)
        ax.legend()

    fig.suptitle("QC Summary", fontsize=14, y=1.02)
    fig.tight_layout()
    _save_figure(fig, save_path)
    return fig
