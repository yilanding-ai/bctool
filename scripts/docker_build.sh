#!/bin/bash
# Build and test the BS-Benchmark Docker image
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SCRIPT_DIR"

echo "========================================"
echo "Building BS-Benchmark Docker Image"
echo "========================================"
echo ""
echo "This will install all 15 aligners via conda."
echo "Estimated time: 15-30 minutes"
echo "Estimated size: 3-5 GB"
echo ""

# Build the image
docker build -t bs-benchmark:latest . 2>&1 | tee docker_build.log

echo ""
echo "========================================"
echo "Build complete! Running verification..."
echo "========================================"

# Verify the build
docker run --rm bs-benchmark:latest --help

echo ""
echo "========================================"
echo "Running demo to verify all tools..."
echo "========================================"

mkdir -p "$SCRIPT_DIR/output"
docker run --rm \
  -v "$SCRIPT_DIR/output:/data/output" \
  bs-benchmark:latest demo -o /data/output/demo_verify

echo ""
echo "========================================"
echo "Docker image ready!"
echo "========================================"
echo ""
echo "Quick commands:"
echo "  # Demo:            docker compose --profile demo run demo"
echo "  # Interactive:     docker compose --profile shell run shell"
echo "  # Run pipeline:    docker compose run bs-benchmark <args>"
echo "  # Manual:          docker run --rm -v \$PWD/output:/data/output bs-benchmark <args>"
echo ""
echo "Example with real data:"
echo "  docker run --rm \\"
echo "    -v \$PWD/input:/data/input:ro \\"
echo "    -v \$PWD/output:/data/output \\"
echo "    -v \$PWD/reference:/data/reference:ro \\"
echo "    bs-benchmark run \\"
echo "      -1 /data/input/sample_R1.fastq.gz \\"
echo "      -2 /data/input/sample_R2.fastq.gz \\"
echo "      -r /data/reference/genome.fa \\"
echo "      -o /data/output/results"
