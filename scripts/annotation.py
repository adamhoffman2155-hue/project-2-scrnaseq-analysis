"""
Cell Type Annotation for Single-Cell RNA-seq Data
Marker-based and reference-based annotation methods
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity


def score_cell_types(adata, marker_dict, score_name="cell_type_score"):
    """
    Score cells based on marker gene expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    marker_dict : dict
        Dictionary mapping cell types to marker genes
        Example: {'T_cell': ['CD3D', 'CD3E'], 'B_cell': ['CD19', 'MS4A1']}
    score_name : str
        Column name for scores

    Returns
    -------
    adata : AnnData
        Updated with cell type scores
    """
    for cell_type, markers in marker_dict.items():
        # Filter markers that exist in data
        markers_present = [m for m in markers if m in adata.var_names]

        if len(markers_present) > 0:
            sc.tl.score_genes(adata, markers_present, score_name=cell_type, use_raw=False)

    return adata


def annotate_by_markers(adata, marker_dict):
    """
    Assign cell types based on highest marker score.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with marker scores
    marker_dict : dict
        Dictionary of cell types and markers

    Returns
    -------
    adata : AnnData
        Updated with cell_type annotation
    """
    # Score cells
    adata = score_cell_types(adata, marker_dict)

    # Get cell type scores
    score_cols = list(marker_dict.keys())
    scores = adata.obs[score_cols]

    # Assign cell type as highest score
    adata.obs["cell_type"] = scores.idxmax(axis=1)
    adata.obs["cell_type_score"] = scores.max(axis=1)

    return adata


def annotate_by_reference(adata, reference_adata, ref_label_key="cell_type"):
    """
    Assign cell types by mapping to reference dataset.

    Parameters
    ----------
    adata : AnnData
        Query dataset
    reference_adata : AnnData
        Reference dataset with known cell types
    ref_label_key : str
        Column in reference.obs with cell type labels

    Returns
    -------
    adata : AnnData
        Query with predicted cell types
    """
    # Find common genes
    common_genes = np.intersect1d(adata.var_names, reference_adata.var_names)

    if len(common_genes) < 100:
        raise ValueError(f"Only {len(common_genes)} common genes found")

    # Subset to common genes
    query_X = adata[:, common_genes].X
    ref_X = reference_adata[:, common_genes].X

    # Compute similarity
    similarity = cosine_similarity(query_X, ref_X)

    # Assign cell type based on nearest neighbor
    nearest_idx = np.argmax(similarity, axis=1)
    adata.obs["cell_type"] = reference_adata.obs[ref_label_key].values[nearest_idx]
    adata.obs["cell_type_confidence"] = np.max(similarity, axis=1)

    return adata


def refine_annotations(adata, cluster_key="leiden", label_key="cell_type"):
    """
    Refine cell type annotations by cluster consensus.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    cluster_key : str
        Column with cluster assignments
    label_key : str
        Column with cell type labels

    Returns
    -------
    adata : AnnData
        Updated with refined annotations
    """
    # Get most common cell type per cluster
    cluster_labels = adata.obs.groupby(cluster_key)[label_key].apply(
        lambda x: x.value_counts().index[0]
    )

    # Update annotations
    adata.obs["cell_type_refined"] = adata.obs[cluster_key].map(cluster_labels)

    return adata


def get_marker_genes(adata, cluster_key="leiden", n_genes=10):
    """
    Get top marker genes for each cell type/cluster.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    cluster_key : str
        Column with cluster/cell type assignments
    n_genes : int
        Number of top genes to return per cluster

    Returns
    -------
    markers_df : pd.DataFrame
        DataFrame with marker genes per cluster
    """
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")

    # Extract results
    markers_list = []
    for cluster in adata.obs[cluster_key].unique():
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster)
        cluster_markers["cluster"] = cluster
        markers_list.append(cluster_markers.head(n_genes))

    markers_df = pd.concat(markers_list, ignore_index=True)

    return markers_df


class CellTypeAnnotator:
    """
    Wrapper class for cell type annotation workflows.
    """

    def __init__(self, marker_dict=None, reference_adata=None):
        """
        Initialize annotator.

        Parameters
        ----------
        marker_dict : dict, optional
            Dictionary of marker genes per cell type
        reference_adata : AnnData, optional
            Reference dataset for mapping
        """
        self.marker_dict = marker_dict
        self.reference_adata = reference_adata

    def annotate(self, adata, method="markers"):
        """
        Annotate cell types.

        Parameters
        ----------
        adata : AnnData
            Query dataset
        method : str
            'markers' or 'reference'

        Returns
        -------
        adata : AnnData
            Annotated dataset
        """
        if method == "markers" and self.marker_dict:
            adata = annotate_by_markers(adata, self.marker_dict)
        elif method == "reference" and self.reference_adata:
            adata = annotate_by_reference(adata, self.reference_adata)
        else:
            raise ValueError(f"Method {method} not available")

        return adata
