"""Tests for scRNA-seq workflow components.

Actually imports and tests the pipeline modules (not generic numpy ops).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Make scripts/ importable
_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(_SCRIPTS))

from generate_synthetic_data import generate_synthetic_data  # noqa: E402
from qc_metrics import (  # noqa: E402
    calculate_qc_metrics,
    detect_outliers,
    filter_cells,
    filter_genes,
)
from annotation import (  # noqa: E402
    score_cell_types,
    annotate_by_markers,
    refine_annotations,
    CellTypeAnnotator,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def synthetic_adata():
    """Small reproducible synthetic dataset."""
    return generate_synthetic_data(n_cells=200, n_genes=500, n_cell_types=3, random_seed=0)


@pytest.fixture
def marker_dict():
    """Marker dictionary using genes that exist in the synthetic data."""
    return {
        "TypeA": [f"Gene_{i}" for i in range(0, 20)],
        "TypeB": [f"Gene_{i}" for i in range(20, 40)],
        "TypeC": [f"Gene_{i}" for i in range(40, 60)],
    }


# ---------------------------------------------------------------------------
# generate_synthetic_data
# ---------------------------------------------------------------------------

def test_synthetic_data_shape(synthetic_adata):
    assert synthetic_adata.n_obs == 200
    assert synthetic_adata.n_vars == 500


def test_synthetic_data_has_qc_columns(synthetic_adata):
    for col in ("n_genes_by_counts", "total_counts", "pct_mt", "pct_rb"):
        assert col in synthetic_adata.obs.columns


def test_synthetic_data_has_cell_types(synthetic_adata):
    assert "cell_type" in synthetic_adata.obs.columns
    assert synthetic_adata.obs["cell_type"].nunique() == 3


def test_synthetic_data_mt_range(synthetic_adata):
    pct_mt = synthetic_adata.obs["pct_mt"].values
    assert (pct_mt >= 0).all() and (pct_mt <= 100).all()


# ---------------------------------------------------------------------------
# qc_metrics
# ---------------------------------------------------------------------------

def test_filter_cells_respects_thresholds(synthetic_adata):
    adata = synthetic_adata.copy()
    n_before = adata.n_obs
    adata = filter_cells(adata, min_genes=0, max_genes=10_000, max_mt=100, max_rb=100)
    assert adata.n_obs == n_before  # permissive filter removes nothing


def test_filter_cells_removes_with_strict_mt(synthetic_adata):
    adata = synthetic_adata.copy()
    adata = filter_cells(adata, min_genes=0, max_genes=10_000, max_mt=2.0, max_rb=100)
    assert adata.n_obs <= synthetic_adata.n_obs  # some cells removed
    assert (adata.obs["pct_mt"] <= 2.0).all()


def test_filter_genes_reduces_var_count(synthetic_adata):
    adata = synthetic_adata.copy()
    n_before = adata.n_vars
    adata = filter_genes(adata, min_cells=5, min_counts=5)
    assert adata.n_vars <= n_before


def test_detect_outliers_returns_bool_mask(synthetic_adata):
    mask = detect_outliers(synthetic_adata, n_mads=3)
    assert mask.dtype == bool
    assert mask.shape == (synthetic_adata.n_obs,)


# ---------------------------------------------------------------------------
# annotation
# ---------------------------------------------------------------------------

def test_score_cell_types_adds_score_columns(synthetic_adata, marker_dict):
    adata = synthetic_adata.copy()
    adata = score_cell_types(adata, marker_dict)
    for ct in marker_dict:
        assert ct in adata.obs.columns


def test_annotate_by_markers_assigns_types(synthetic_adata, marker_dict):
    adata = synthetic_adata.copy()
    adata = annotate_by_markers(adata, marker_dict)
    assert "cell_type" in adata.obs.columns
    assigned = set(adata.obs["cell_type"].unique())
    assert assigned.issubset(set(marker_dict.keys()))


def test_refine_annotations_requires_cluster_column(synthetic_adata, marker_dict):
    adata = synthetic_adata.copy()
    adata = annotate_by_markers(adata, marker_dict)
    # Fake cluster assignments
    adata.obs["leiden"] = (np.arange(adata.n_obs) % 4).astype(str)
    adata = refine_annotations(adata, cluster_key="leiden", label_key="cell_type")
    assert "cell_type_refined" in adata.obs.columns
    assert len(adata.obs["cell_type_refined"]) == adata.n_obs


def test_cell_type_annotator_class(synthetic_adata, marker_dict):
    adata = synthetic_adata.copy()
    annotator = CellTypeAnnotator(marker_dict=marker_dict)
    adata = annotator.annotate(adata, method="markers")
    assert "cell_type" in adata.obs.columns


def test_cell_type_annotator_rejects_unknown_method(synthetic_adata, marker_dict):
    annotator = CellTypeAnnotator(marker_dict=marker_dict)
    with pytest.raises(ValueError):
        annotator.annotate(synthetic_adata, method="nonsense")
