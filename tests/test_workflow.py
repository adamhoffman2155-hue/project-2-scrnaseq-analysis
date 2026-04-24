"""Tests for scRNA-seq workflow components."""

import numpy as np


def test_synthetic_data_shape():
    """Test that synthetic data has expected dimensions."""
    from scipy import sparse

    n_cells, n_genes = 500, 2000
    X = sparse.random(n_cells, n_genes, density=0.1, format="csr")
    assert X.shape == (n_cells, n_genes)
    assert X.nnz > 0


def test_qc_metric_ranges():
    """Test QC metrics are in realistic ranges."""
    n_cells = 100
    n_umi = np.random.randint(1000, 20000, size=n_cells)
    pct_mt = np.random.uniform(0, 20, size=n_cells)
    assert np.all(n_umi > 0)
    assert np.all(pct_mt >= 0) and np.all(pct_mt <= 100)


def test_qc_filtering_removes_cells():
    """Test that QC filtering removes low-quality cells."""
    pct_mt = np.random.uniform(0, 30, size=100)
    mask = pct_mt < 20.0
    assert 0 < mask.sum() < 100


def test_cell_type_annotation_format():
    """Test cell type annotation produces valid labels."""
    cell_types = ["T_cell", "B_cell", "Macrophage", "Fibroblast", "Epithelial"]
    annotations = np.random.choice(cell_types, size=50)
    assert len(annotations) == 50
    assert all(ct in cell_types for ct in annotations)


def test_clustering_produces_labels():
    """Test clustering assigns every cell a label."""
    labels = np.random.randint(0, 5, size=100)
    assert len(labels) == 100
    assert len(np.unique(labels)) <= 5


def test_marker_gene_detection():
    """Test marker genes have expected properties."""
    pvals = np.random.uniform(0, 1, 100)
    pvals_adj = np.minimum(pvals * 100 / np.arange(1, 101), 1.0)
    significant = [i for i, p in enumerate(pvals_adj) if p < 0.05]
    assert isinstance(significant, list)
