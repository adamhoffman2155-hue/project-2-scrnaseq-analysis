"""
Generate Synthetic scRNA-seq Data
Creates a realistic sparse count matrix with known cell types and QC metrics.
"""

import numpy as np
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix


def generate_synthetic_data(
    n_cells=500,
    n_genes=2000,
    n_cell_types=5,
    random_seed=42,
    output_path=None,
):
    """
    Create a synthetic scRNA-seq AnnData object.

    Parameters
    ----------
    n_cells : int
        Number of cells to simulate.
    n_genes : int
        Number of genes.
    n_cell_types : int
        Number of distinct cell types.
    random_seed : int
        Random seed for reproducibility.
    output_path : str or None
        If provided, save as h5ad.

    Returns
    -------
    adata : AnnData
    """
    rng = np.random.RandomState(random_seed)

    # --- cell type assignments ---
    cell_type_names = [f"CellType_{i}" for i in range(n_cell_types)]
    labels = rng.choice(cell_type_names, size=n_cells)

    # --- base expression (sparse Poisson) ---
    lam_base = 0.3  # low background rate
    counts = rng.poisson(lam=lam_base, size=(n_cells, n_genes)).astype(np.float32)

    # --- marker genes: each cell type gets a block of up-regulated genes ---
    markers_per_type = 20
    marker_gene_map = {}
    for idx, ct in enumerate(cell_type_names):
        start = idx * markers_per_type
        end = start + markers_per_type
        marker_gene_map[ct] = list(range(start, end))
        mask = labels == ct
        counts[mask, start:end] += rng.poisson(lam=5, size=(mask.sum(), markers_per_type))

    # --- realistic total UMI (1000-20000) ---
    target_umi = rng.randint(1000, 20001, size=n_cells).astype(np.float32)
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    counts = (counts / row_sums * target_umi[:, None]).astype(np.float32)
    counts = np.round(counts).astype(np.float32)

    # --- gene names ---
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    # Add a few "mitochondrial" genes
    for i in range(n_genes - 5, n_genes):
        gene_names[i] = f"MT-G{i}"

    # --- build AnnData ---
    sparse_counts = csr_matrix(counts)
    adata = ad.AnnData(X=sparse_counts)
    adata.var_names = gene_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]
    adata.obs["cell_type"] = labels

    # --- QC-like metadata ---
    adata.obs["n_genes_by_counts"] = np.array((counts > 0).sum(axis=1)).flatten()
    adata.obs["total_counts"] = np.array(counts.sum(axis=1)).flatten()

    # pct_mt: fraction from last 5 genes (MT-*)
    mt_counts = counts[:, -5:].sum(axis=1)
    total = counts.sum(axis=1)
    total[total == 0] = 1
    adata.obs["pct_mt"] = (mt_counts / total * 100).astype(np.float32)
    # Clamp pct_mt to 1-15% range for realism
    adata.obs["pct_mt"] = adata.obs["pct_mt"].clip(1.0, 15.0)

    adata.obs["pct_rb"] = rng.uniform(5, 30, size=n_cells).astype(np.float32)

    # Store marker map in uns
    adata.uns["marker_gene_map"] = {
        ct: [gene_names[g] for g in idxs] for ct, idxs in marker_gene_map.items()
    }

    if output_path is not None:
        adata.write_h5ad(output_path)
        print(f"Saved synthetic data to {output_path}")

    return adata


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate synthetic scRNA-seq data")
    parser.add_argument("-o", "--output", default="synthetic_data.h5ad",
                        help="Output h5ad path")
    parser.add_argument("--n-cells", type=int, default=500)
    parser.add_argument("--n-genes", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    adata = generate_synthetic_data(
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        random_seed=args.seed,
        output_path=args.output,
    )
    print(f"Generated: {adata.n_obs} cells x {adata.n_vars} genes")
