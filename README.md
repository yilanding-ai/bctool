# BCTool - Bisulfite Conversion Alignment & Benchmarking Toolkit

Multi-algorithm methylation alignment benchmark with a **built-in conversion-aware aligner (ConvSeed)** supporting all 12 single-base substitution types.

## Features

- **Built-in ConvSeed aligner**: 2-base seed + 2-bit bitmask + banded Smith-Waterman, covers all 12 conversion types without external tools
- **15 external alignment tools**: Bwa-meth, BSBolt, BSMAP, Walt, Abismal, Batmeth2, Basal, HISAT-3n, HISAT-3n (repeat), Bismark-bwt2-e2e, Bismark-his2, BSseeker2-bwt/soap2/bwt2-e2e/bwt2-local
- **12 conversion types**: C>T, T>C, A>G, G>A, A>C, C>A, G>T, T>G, A>T, T>A, C>G, G>C
- **Smart error correction**: 7 strategies across 2 phases
  - SAM-level: mq, clip, pair, unconverted, model (ML-based)
  - BED-level: consensus (hard-threshold), context (motif-aware Bayesian)
- **Motif detection**: Auto-discovers methylation-enriched sequence motifs from reference
- **ML-based error detection**: Logistic regression on 9 SAM features, pure numpy
- **Algorithm selector**: Auto-picks correction params by read length (4 tiers)
- **Simulated data**: Built-in simulator with known ground truth
- **Real data support**: Extract subsets from real FASTQ files, adapter trimming, UMI processing, QC
- **Single-end & Paired-end**: Both modes supported
- **Multi-sample batch analysis**: Compare aligner performance across samples
- **Web UI**: Interactive Streamlit interface with Plotly charts
- **Docker**: Ready-to-deploy container with all tools

## Installation

```bash
# From source
git clone https://github.com/yilanding-ai/bctool.git
cd bctool
pip install -r requirements.txt
pip install .

# Verify
bctool --help
```

## Quick Start

```bash
# Demo with built-in ConvSeed aligner (no external tools)
bctool run --simulate --aligners conversion_aware -o ./results

# Run all built-in mock aligners + ConvSeed
bctool run --simulate -o ./results

# Demo with error correction
bctool run --simulate --error-correct --correction-strategies mq,clip,context -o ./results

# Web UI
bctool web
```

## ConvSeed: Built-in Conversion-Aware Aligner

ConvSeed is a purpose-built aligner that combines 3 complementary techniques:

| Component | Technique | Detail |
|---|---|---|
| **Seed encoding** | 3-base conversion-aware | 4 bases → 3 states: conversion pair shares one state, others distinct; k=8-16 auto-adapted by read length |
| **Seed lookup** | Rolling hash (base-3) | O(1) hash table, 3^k entries eliminates collision for genome ≤ 1Mbp |
| **Alignment** | Banded SW (vectorized) | Affine gap (open=-4, extend=-1), match=+1, conversion-match=0, mismatch=-1; numpy-vectorized inner loop |

Registration name: `conversion_aware`. Automatically handles all 12 conversion types with zero configuration.

## Error Correction

| Phase | Strategy | Description | Key Parameter |
|---|---|---|---|---|
| SAM | `mq` | Remove low-MAPQ reads | `--correction-min-mq` (default: 20) |
| SAM | `clip` | Remove excessive soft-clip | `max_clip_pct: 50` |
| SAM | `pair` | Remove discordant PE reads | TLEN > 1000bp or wrong chr |
| SAM | `unconverted` | Remove reads with too many unconverted target bases; writes rejected reads to `.fastq` | `--correction-max-unconverted` (default: 3) |
| SAM | `model` | Logistic regression on 9 features; trained against ground truth | `--correction-model-threshold` (default: 0.5) |
| BED | `consensus` | Hard threshold: ratio ≥0.7→1.0, ≤0.3→0.0 | All read lengths |
| BED | `context` | Beta-Binomial with motif-specific Bayesian priors; extracts ±N flank from reference | <80bp: simple<br>80-149: flank=1<br>150-499: flank=2<br>≥500: flank=3 |

## CLI

```text
bctool [OPTIONS] COMMAND [ARGS]...

Commands:
  run              Run the full benchmarking pipeline
  simulate         Generate simulated bisulfite/converted sequencing data
  extract          Extract a subset of reads from FASTQ files
  trim             Trim adapters from FASTQ files
  umi              Detect and process UMIs
  qc               Quality control analysis
  report           Generate report from completed benchmark results
  batch            Multi-sample batch analysis
  list-aligners    List all supported aligners
  web              Launch Streamlit web interface
```

## Examples

```bash
# Run ConvSeed on C>T conversion
bctool run --simulate --conversion ct --aligners conversion_aware -o ./ct_results

# All 12 conversions in one batch
for conv in ct tc ag ga ac ca gt tg at ta cg gc; do
  bctool run --simulate --conversion $conv --aligners conversion_aware -o ./${conv}_results
done

# Full pipeline with error correction
bctool run --simulate --error-correct \
  --correction-strategies mq,clip,unconverted,context \
  --correction-min-mq 30 -o ./corrected_results

# Batch demo (4 samples, 2 conversion types)
bctool batch-demo

# Real data
bctool run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -r reference.fa --conversion ct -o ./results --threads 16
```

## Conversion Types

| Tag | Conversion | Description |
|---|---|---|
| ct | C→T | Bisulfite (WGBS) |
| tc | T→C | Reverse bisulfite |
| ag | A→G | A-to-G conversion |
| ga | G→A | G-to-A conversion |
| ac | A→C | A-to-C conversion |
| ca | C→A | C-to-A conversion |
| gt | G→T | G-to-T conversion |
| tg | T→G | T-to-G conversion |
| at | A→T | A-to-T conversion |
| ta | T→A | T-to-A conversion |
| cg | C→G | C-to-G conversion |
| gc | G→C | G-to-C conversion |

## Output Structure

```
results/
+-- comparison_results.json       # Full results in JSON
+-- comparison_results.csv        # Results in CSV
+-- report/
|   +-- report.html               # Interactive HTML report
|   +-- summary.txt               # Text summary
|   +-- f1_score.png              # F1 score bar chart
|   +-- precision.png             # Precision bar chart
|   +-- recall.png                # Recall bar chart
|   +-- accuracy.png              # Accuracy bar chart
|   +-- heatmap.png               # Performance heatmap
+-- qc/                           # QC reports (if enabled)
+-- work/
    +-- simulated/                # Simulated data (if --simulate)
    +-- conversion_aware/         # Per-aligner results
```

## License

MIT
