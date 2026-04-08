FROM continuumio/miniconda3:23.5.2-0

LABEL maintainer="Adam Hoffman <adamhoffman21@hotmail.ca>"
LABEL description="Single-Cell RNA-seq Analysis Workflow (Seurat + Scanpy)"

WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file
COPY environment.yml .

# Create conda environment
RUN conda env create -f environment.yml && \
    conda clean --all --yes

# Activate environment
ENV PATH /opt/conda/envs/scrnaseq-analysis/bin:$PATH
SHELL ["/bin/bash", "-c"]

# Verify installations
RUN conda run -n scrnaseq-analysis \
    python -c "import scanpy; print(f'Scanpy {scanpy.__version__}')" && \
    Rscript -e "library(Seurat); cat('Seurat loaded successfully\n')"

# Copy project files
COPY scripts/ scripts/
COPY config/ config/
COPY tests/ tests/

# Create data directories
RUN mkdir -p data/{raw,processed,metadata} results/{qc,figures,tables,objects}

# Set entrypoint
ENTRYPOINT ["python"]
CMD ["scripts/run_workflow.py"]
