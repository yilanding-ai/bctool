#!/bin/bash
# Install all 15 alignment tools via conda
# Usage: bash install_tools.sh

set -e

echo "========================================"
echo "Installing 14 WGBS alignment tools"
echo "========================================"

# Create conda environment
ENV_NAME="bs-benchmark"
if ! conda env list | grep -q "^$ENV_NAME "; then
    echo "Creating conda environment: $ENV_NAME"
    conda create -y -n $ENV_NAME python=3.9
fi

eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Install Python dependencies
pip install click pyyaml pandas numpy biopython matplotlib seaborn rich

# Install core aligners
echo "Installing core aligners..."
conda install -y -c bioconda \
    bwameth \
    bsbolt \
    bsmap \
    bismark \
    hisat2 hisat-3n \
    bs-seeker2 \
    samtools \
    bedtools

# Install additional aligners
echo "Installing additional aligners..."

# Walt
pip install walt-bio

# Abismal
conda install -y -c bioconda abismal

# Batmeth2
pip install batmeth2

# BSseeker2 components
pip install bs-seeker2

echo ""
echo "Installation complete!"
echo "Activate environment: conda activate $ENV_NAME"
echo ""
echo "Verify installation: bs-benchmark run --dry-run -1 test.fq"

# Check availability
echo ""
echo "=== Tool Availability Check ==="
for tool in bwameth.py bsbolt bsmap walt abismal batmeth2 hisat-3n bismark bs_seeker2-align.py; do
    if command -v $tool &> /dev/null; then
        echo "[OK] $tool"
    else
        echo "[MISSING] $tool"
    fi
done
