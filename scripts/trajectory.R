#!/usr/bin/env Rscript

# Monocle3 Trajectory Analysis for scRNA-seq Data
# Converts Seurat object to CDS, learns trajectory, and orders cells in pseudotime.

library(Seurat)
library(monocle3)
library(ggplot2)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript trajectory.R <input_seurat_rds> <output_dir> [root_cell_type]\n")
  cat("  input_seurat_rds  Path to a saved Seurat RDS object\n")
  cat("  output_dir        Directory for output files\n")
  cat("  root_cell_type    (optional) Cell type to use as trajectory root\n")
  quit(status = 1)
}

input_path   <- args[1]
output_dir   <- args[2]
root_type    <- if (length(args) >= 3) args[3] else NULL

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helper: convert Seurat -> cell_data_set
# ---------------------------------------------------------------------------
seurat_to_cds <- function(seurat_obj) {
  # Extract data
  expression_matrix <- GetAssayData(seurat_obj, slot = "counts")
  cell_metadata     <- seurat_obj@meta.data
  gene_metadata     <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names       = rownames(expression_matrix)
  )

  cds <- new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata   = cell_metadata,
    gene_metadata   = gene_metadata
  )
  return(cds)
}

# ---------------------------------------------------------------------------
# Helper: pick root cells
# ---------------------------------------------------------------------------
get_root_cells <- function(cds, cell_type_col = "cell_type", root_type = NULL) {
  if (!is.null(root_type) && cell_type_col %in% colnames(colData(cds))) {
    root_cells <- colnames(cds)[colData(cds)[[cell_type_col]] == root_type]
    if (length(root_cells) == 0) {
      warning("No cells found for root type '", root_type, "'. Using first principal node.")
      return(NULL)
    }
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes  <- igraph::V(principal_graph(cds)[["UMAP"]])$name
    root_pr_node   <- names(which.max(table(
      closest_vertex[root_cells, ]
    )))
    return(root_pr_node)
  }
  return(NULL)
}

# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
cat("=== Monocle3 Trajectory Analysis ===\n")

cat("Loading Seurat object:", input_path, "\n")
seurat_obj <- readRDS(input_path)

cat("Converting to cell_data_set...\n")
cds <- seurat_to_cds(seurat_obj)

cat("Preprocessing (normalization + PCA)...\n")
cds <- preprocess_cds(cds, num_dim = 30, method = "PCA")

cat("Reducing dimensions (UMAP)...\n")
cds <- reduce_dimension(cds, reduction_method = "UMAP")

cat("Clustering cells...\n")
cds <- cluster_cells(cds, reduction_method = "UMAP")

cat("Learning trajectory graph...\n")
cds <- learn_graph(cds, use_partition = TRUE)

cat("Ordering cells in pseudotime...\n")
root_node <- get_root_cells(cds, root_type = root_type)
if (!is.null(root_node)) {
  cds <- order_cells(cds, root_pr_nodes = root_node)
} else {
  cds <- order_cells(cds)
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
cat("Generating trajectory plots...\n")

# Trajectory colored by cell type
if ("cell_type" %in% colnames(colData(cds))) {
  p1 <- plot_cells(cds, color_cells_by = "cell_type",
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
    ggtitle("Trajectory colored by cell type")
  ggsave(file.path(output_dir, "trajectory_celltype.png"), p1,
         width = 8, height = 6, dpi = 300)
}

# Trajectory colored by pseudotime
p2 <- plot_cells(cds, color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 label_leaves = FALSE, label_branch_points = FALSE) +
  ggtitle("Trajectory colored by pseudotime")
ggsave(file.path(output_dir, "trajectory_pseudotime.png"), p2,
       width = 8, height = 6, dpi = 300)

# Trajectory colored by cluster partition
p3 <- plot_cells(cds, color_cells_by = "partition",
                 label_groups_by_cluster = FALSE) +
  ggtitle("Trajectory colored by partition")
ggsave(file.path(output_dir, "trajectory_partition.png"), p3,
       width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
cat("Saving outputs...\n")
saveRDS(cds, file.path(output_dir, "monocle3_cds.rds"))

# Export pseudotime values
pseudotime_df <- data.frame(
  cell     = colnames(cds),
  pseudotime = pseudotime(cds)
)
write.csv(pseudotime_df, file.path(output_dir, "pseudotime_values.csv"),
          row.names = FALSE)

cat("Done. Results saved to:", output_dir, "\n")
