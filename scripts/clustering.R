#!/usr/bin/env Rscript

# Seurat-based Clustering and Dimensionality Reduction

library(Seurat)
library(tidyverse)
library(ggplot2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Create Seurat object from count matrix
create_seurat_object <- function(counts_path, project_name = "scRNA") {
  cat("Loading count matrix...\n")
  
  # Read counts (assuming CSV or MTX format)
  if (grepl("\\.csv", counts_path)) {
    counts <- read.csv(counts_path, row.names = 1)
  } else if (grepl("\\.mtx", counts_path)) {
    counts <- Matrix::readMM(counts_path)
  } else {
    stop("Unsupported file format")
  }
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = 3,
    min.features = 200
  )
  
  return(seurat_obj)
}

# Quality control and filtering
perform_qc <- function(seurat_obj, max_mt = 20) {
  cat("Calculating QC metrics...\n")
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < max_mt
  )
  
  cat(paste("After QC filtering:", ncol(seurat_obj), "cells\n"))
  
  return(seurat_obj)
}

# Normalization and scaling
normalize_and_scale <- function(seurat_obj) {
  cat("Normalizing and scaling...\n")
  
  # Normalize
  seurat_obj <- NormalizeData(seurat_obj)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  
  # Scale
  seurat_obj <- ScaleData(seurat_obj)
  
  return(seurat_obj)
}

# Dimensionality reduction
reduce_dimensions <- function(seurat_obj, n_pcs = 50) {
  cat("Performing PCA...\n")
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs)
  
  cat("Computing UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  
  cat("Computing t-SNE...\n")
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
  
  return(seurat_obj)
}

# Clustering
perform_clustering <- function(seurat_obj, resolution = 0.5) {
  cat("Finding neighbors...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  
  cat(paste("Clustering with resolution:", resolution, "\n"))
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  return(seurat_obj)
}

# Find cluster markers
find_markers <- function(seurat_obj) {
  cat("Finding cluster markers...\n")
  
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  return(markers)
}

# Main workflow
main <- function(input_path, output_dir, resolution = 0.5) {
  cat("=== Seurat Clustering Workflow ===\n")

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load Seurat object (RDS or count matrix)
  cat("Loading data from:", input_path, "\n")
  if (grepl("\\.rds$", input_path, ignore.case = TRUE)) {
    seurat_obj <- readRDS(input_path)
  } else {
    seurat_obj <- create_seurat_object(input_path)
  }

  # QC filtering
  seurat_obj <- perform_qc(seurat_obj)

  # Find variable features
  cat("Finding variable features...\n")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

  # Scale data
  cat("Scaling data...\n")
  seurat_obj <- ScaleData(seurat_obj)

  # PCA
  cat("Running PCA...\n")
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)

  # Neighbor graph
  cat("Finding neighbors...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

  # Clustering
  cat(paste("Clustering with resolution:", resolution, "\n"))
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

  # UMAP
  cat("Running UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

  # Save clustered Seurat object
  rds_path <- file.path(output_dir, "clustered_seurat.rds")
  saveRDS(seurat_obj, rds_path)
  cat("Saved clustered object to:", rds_path, "\n")

  # Save UMAP plot
  umap_path <- file.path(output_dir, "umap_clusters.png")
  p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle("UMAP - Seurat Clusters")
  ggsave(umap_path, plot = p, width = 8, height = 6, dpi = 300)
  cat("Saved UMAP plot to:", umap_path, "\n")

  cat("Clustering workflow complete.\n")
  return(seurat_obj)
}

# Run when executed as a script
if (!interactive()) {
  if (length(args) < 2) {
    cat("Usage: Rscript clustering.R <input_path> <output_dir> [resolution]\n")
    cat("  input_path   Path to Seurat RDS or count matrix (CSV/MTX)\n")
    cat("  output_dir   Directory for output files\n")
    cat("  resolution   (optional) Clustering resolution, default 0.5\n")
    quit(status = 1)
  }

  input_path <- args[1]
  output_dir <- args[2]
  resolution <- if (length(args) >= 3) as.numeric(args[3]) else 0.5

  main(input_path, output_dir, resolution)
}
