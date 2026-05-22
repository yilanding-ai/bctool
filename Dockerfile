# =============================================================================
# bctool Docker Image
# All 14 WGBS alignment tools + Python pipeline in one container
# =============================================================================

FROM ubuntu:22.04 AS base

LABEL maintainer="bctool" \
      description="bctool: multi-algorithm methylation alignment benchmark"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl ca-certificates \
    build-essential gcc g++ make cmake \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libncurses5-dev libncursesw5-dev \
    python3 python3-pip python3-dev \
    git automake autoconf \
    libhts-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# Stage: Install bioinformatics tools via conda
# =============================================================================
FROM base AS conda-env

# Install Miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
    && bash /miniconda.sh -b -p /opt/conda \
    && rm /miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

# Install all 15 aligners + utilities via conda/bioconda
RUN conda create -y -n bctool python=3.10 && \
    conda install -y -n bctool -c bioconda -c conda-forge \
    # Core aligners from the paper
    bwameth \
    bsbolt \
    bsmap \
    bismark \
    hisat2 \
    hisat-3n \
    bs-seeker2 \
    abismal \
    walt-bio \
    batmeth2 \
    # Utilities
    samtools \
    bedtools \
    cutadapt \
    fastqc \
    multiqc \
    seqtk \
    # Python packages
    && conda clean -afy

# =============================================================================
# Stage: Build final image
# =============================================================================
FROM conda-env AS final

ENV PATH=/opt/conda/envs/bctool/bin:/opt/conda/bin:$PATH
SHELL ["conda", "run", "-n", "bctool", "/bin/bash", "-c"]

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Copy bctool package
COPY . /opt/bctool
WORKDIR /opt/bctool
RUN pip install -e .

# Create data directories
RUN mkdir -p /data/input /data/output /data/reference

# Verify tool installation
RUN echo "=== Tool Verification ===" && \
    for tool in bwameth.py bsbolt bsmap bismark hisat-3n abismal cutadapt seqtk; do \
        which $tool && echo "[OK] $tool" || echo "[MISSING] $tool"; \
    done && \
    python -c "import bctool; print('bctool:', bctool.__version__)"

# Entry point
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bctool", "bctool"]
CMD ["--help"]

# =============================================================================
# Usage:
#   docker build -t bctool .
#
#   # Demo mode (no input files needed)
#   docker run --rm -v $(pwd)/output:/data/output bctool demo -o /data/output
#
#   # Run with real data
#   docker run --rm \
#     -v $(pwd)/input:/data/input \
#     -v $(pwd)/output:/data/output \
#     -v $(pwd)/reference:/data/reference \
#     bctool run \
#       -1 /data/input/sample_R1.fastq.gz \
#       -2 /data/input/sample_R2.fastq.gz \
#       -r /data/reference/genome.fa \
#       -o /data/output
#
#   # Batch mode
#   docker run --rm \
#     -v $(pwd)/batch_config.yaml:/data/config.yaml \
#     -v $(pwd)/output:/data/output \
#     bctool batch run -c /data/config.yaml -o /data/output
#
#   # QC only
#   docker run --rm \
#     -v $(pwd)/input:/data/input \
#     -v $(pwd)/output:/data/output \
#     bctool qc /data/input/sample_R1.fastq -o /data/output
# =============================================================================
