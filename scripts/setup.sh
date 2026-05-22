#!/bin/bash
# One-command setup for BS-Benchmark
# Usage: bash setup.sh [--docker|--conda]

set -e

echo "========================================"
echo "BS-Benchmark Setup"
echo "========================================"

MODE="${1:---conda}"

case "$MODE" in
    --docker)
        echo "Mode: Docker container"
        echo ""
        echo "Building Docker image (15-30 min, ~4GB)..."
        bash "$(dirname "$0")/docker_build.sh"
        ;;

    --conda)
        echo "Mode: Local conda environment"
        echo ""

        # Check conda
        if ! command -v conda &> /dev/null; then
            echo "ERROR: conda not found."
            echo "Install Miniconda first:"
            echo "  https://docs.conda.io/en/latest/miniconda.html"
            exit 1
        fi

        ENV_NAME="${CONDA_ENV:-bs-benchmark}"

        echo "Creating conda environment: $ENV_NAME"
        conda create -y -n $ENV_NAME python=3.10

        # Install tools
        echo "Installing 15 alignment tools..."
        conda install -y -n $ENV_NAME -c bioconda -c conda-forge \
            bwameth \
            bsbolt \
            bsmap \
            bismark \
            hisat2 \
            hisat-3n \
            bs-seeker2 \
            abismal \
            samtools \
            bedtools \
            cutadapt \
            fastqc \
            seqtk

        # Pip-based tools
        echo "Installing pip-based tools..."
        conda run -n $ENV_NAME pip install walt-bio batmeth2

        # Python packages
        echo "Installing Python dependencies..."
        conda run -n $ENV_NAME pip install -r "$(dirname "$0")/../requirements.txt"

        # Install bs-benchmark package
        echo "Installing BS-Benchmark..."
        conda run -n $ENV_NAME pip install -e "$(dirname "$0")/.."

        echo ""
        echo "========================================"
        echo "Setup complete!"
        echo "========================================"
        echo "Activate:  conda activate $ENV_NAME"
        echo "Quick test: bs-benchmark demo -o ./demo_results"
        echo "Dry run:   bs-benchmark run --dry-run -1 test.fq"
        ;;

    *)
        echo "Usage: bash setup.sh [--docker|--conda]"
        exit 1
        ;;
esac
