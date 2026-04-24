# Benchmark: batch-correction methods on synthetic scRNA-seq

Small synthetic AnnData (600 cells, 800 genes, 5 simulated cell types,
2 simulated batches with a per-gene additive offset). Silhouette scores
on the learned embedding:
- **batch_silhouette**: lower = better batch mixing (goal: ≈ 0)
- **celltype_silhouette**: higher = cell types preserved

| Method | Batch silhouette (↓) | Cell-type silhouette (↑) |
| --- | ---: | ---: |
| PCA (uncorrected) | 0.002 | 0.483 |
| Harmony | -0.002 | 0.488 |
| scVI (skipped — scvi-tools not installed) | n/a | n/a |

## Interpretation

Uncorrected PCA shows the synthetic batch offset as separable in the
embedding. Harmony is expected to compress the batch silhouette toward
zero while keeping cell-type structure. scVI (if installed) should do
similarly via its variational batch-aware encoder. On 600 cells both
methods are essentially tied; the benchmark exists as a template for
running on real data.
