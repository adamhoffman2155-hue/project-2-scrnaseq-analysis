#!/usr/bin/env Rscript

"""
Seurat-based Clustering and Dimensionality Reduction
"""

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
    counts <- readMM(counts_path)
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
main <- function() {
  # Example usage
  cat("=== Seurat Clustering Workflow ===\n")
  
  # This would be called with actual data
  # seurat_obj <- create_seurat_object("data/raw/counts.csv")
  # seurat_obj <- perform_qc(seurat_obj)
  # seurat_obj <- normalize_and_scale(seurat_obj)
  # seurat_obj <- reduce_dimensions(seurat_obj)
  # seurat_obj <- perform_clustering(seurat_obj, resolution = 0.5)
  # markers <- find_markers(seurat_obj)
  
  cat("Workflow functions defined. Use in interactive R session or Jupyter.\n")
}

# Export functions for use in notebooks
if (!interactive()) {
  main()
}
