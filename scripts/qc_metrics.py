"""
Quality Control Metrics for Single-Cell RNA-seq Data
Calculate and visualize QC metrics
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import median_abs_deviation


def calculate_qc_metrics(adata, mt_pattern="^MT-", rb_pattern="^RPS|^RPL"):
    """
    Calculate QC metrics for single-cell data.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    mt_pattern : str
        Regex pattern for mitochondrial genes
    rb_pattern : str
        Regex pattern for ribosomal genes
    
    Returns
    -------
    adata : AnnData
        Updated AnnData with QC metrics in .obs
    """
    import scanpy as sc
    
    # Calculate standard metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    # Calculate mitochondrial percentage
    mt_genes = adata.var_names.str.contains(mt_pattern, case=False)
    adata.obs['pct_mt'] = np.sum(adata[:, mt_genes].X, axis=1) / np.sum(adata.X, axis=1) * 100
    
    # Calculate ribosomal percentage
    rb_genes = adata.var_names.str.contains(rb_pattern, case=False)
    adata.obs['pct_rb'] = np.sum(adata[:, rb_genes].X, axis=1) / np.sum(adata.X, axis=1) * 100
    
    return adata


def detect_outliers(adata, n_mads=3):
    """
    Detect QC outliers using MAD (Median Absolute Deviation).
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics
    n_mads : int
        Number of MADs for outlier threshold
    
    Returns
    -------
    outlier_mask : np.ndarray
        Boolean array indicating outliers
    """
    qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_mt']
    
    outlier_mask = np.zeros(adata.n_obs, dtype=bool)
    
    for metric in qc_metrics:
        if metric in adata.obs.columns:
            values = adata.obs[metric].values
            median = np.median(values)
            mad = median_abs_deviation(values)
            threshold_lower = median - n_mads * mad
            threshold_upper = median + n_mads * mad
            
            outliers = (values < threshold_lower) | (values > threshold_upper)
            outlier_mask |= outliers
    
    return outlier_mask


def plot_qc_metrics(adata, save_path=None):
    """
    Generate QC diagnostic plots.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with QC metrics
    save_path : str, optional
        Path to save figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # nGenes vs nCounts
    axes[0, 0].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], 
                       alpha=0.5, s=10)
    axes[0, 0].set_xlabel('Total Counts')
    axes[0, 0].set_ylabel('Number of Genes')
    axes[0, 0].set_title('Gene Count vs Total Counts')
    
    # Mitochondrial percentage
    axes[0, 1].hist(adata.obs['pct_mt'], bins=50, edgecolor='black')
    axes[0, 1].set_xlabel('% Mitochondrial Genes')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Mitochondrial Gene Distribution')
    
    # Ribosomal percentage
    axes[1, 0].hist(adata.obs['pct_rb'], bins=50, edgecolor='black')
    axes[1, 0].set_xlabel('% Ribosomal Genes')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Ribosomal Gene Distribution')
    
    # Total counts distribution
    axes[1, 1].hist(np.log10(adata.obs['total_counts'] + 1), bins=50, edgecolor='black')
    axes[1, 1].set_xlabel('log10(Total Counts)')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Total Counts Distribution')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def filter_cells(adata, min_genes=200, max_genes=5000, max_mt=20, max_rb=50):
    """
    Filter low-quality cells.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    min_genes : int
        Minimum number of genes per cell
    max_genes : int
        Maximum number of genes per cell
    max_mt : float
        Maximum mitochondrial percentage
    max_rb : float
        Maximum ribosomal percentage
    
    Returns
    -------
    adata_filtered : AnnData
        Filtered AnnData object
    """
    initial_n = adata.n_obs
    
    adata = adata[
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['pct_mt'] <= max_mt) &
        (adata.obs['pct_rb'] <= max_rb)
    ]
    
    print(f"Filtered {initial_n - adata.n_obs} cells ({100*(initial_n-adata.n_obs)/initial_n:.1f}%)")
    print(f"Remaining: {adata.n_obs} cells")
    
    return adata


def filter_genes(adata, min_cells=3, min_counts=3):
    """
    Filter lowly-expressed genes.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    min_cells : int
        Minimum number of cells expressing gene
    min_counts : int
        Minimum total counts for gene
    
    Returns
    -------
    adata_filtered : AnnData
        Filtered AnnData object
    """
    initial_n = adata.n_vars
    
    # Count cells expressing each gene
    expressed = np.sum(adata.X > 0, axis=0).A1 if hasattr(adata.X, 'A1') else np.sum(adata.X > 0, axis=0)
    total_counts = np.sum(adata.X, axis=0).A1 if hasattr(adata.X, 'A1') else np.sum(adata.X, axis=0)
    
    adata = adata[:, (expressed >= min_cells) & (total_counts >= min_counts)]
    
    print(f"Filtered {initial_n - adata.n_vars} genes ({100*(initial_n-adata.n_vars)/initial_n:.1f}%)")
    print(f"Remaining: {adata.n_vars} genes")
    
    return adata
