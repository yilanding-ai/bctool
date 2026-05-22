"""Error correction module for post-alignment read filtering and correction.

Strategies:
  - mq:   Remove reads below mapping quality threshold
  - clip: Remove/examine reads with excessive soft-clipping
  - pair: Filter discordant paired-end reads (wrong orientation/distance)
  - unconverted: Remove reads with too many unconverted target bases
  - consensus: Correct methylation calls based on per-site consensus
  - context: Motif-aware Bayesian correction (auto-selects algorithm by read length)
"""

import re
from pathlib import Path
from collections import defaultdict


# IUPAC ambiguity â†?regex character class
_IUPAC_RE = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def _iupac_to_regex(pattern):
    """Convert IUPAC motif pattern to regex string."""
    return "".join(_IUPAC_RE.get(c.upper(), c) for c in pattern)


class AlgorithmSelector:
    """Select correction algorithm based on read length and data characteristics.

    Read length tiers:
      < 80 bp     â†?simple       (too few target sites, minimal context)
      80â€?49 bp   â†?context      (light: flank=1, prior=5)
      150â€?99 bp  â†?context      (standard: flank=2, prior=10)
      >= 500 bp   â†?haplotype    (read-level HMM; falls back to context if sparse)
    """

    @staticmethod
    def select(read_length, target_density=None):
        """Return (algorithm_name, recommended_flank, recommended_prior_strength)."""
        if read_length < 80:
            return "simple", 0, 5
        if read_length < 150:
            return "context", 1, 5
        if read_length < 500:
            return "context", 2, 10
        if target_density is not None and target_density < 0.05:
            return "context", 2, 10  # too sparse for haplotype
        return "haplotype", 3, 10


class MotifDetector:
    """Auto-detect methylation-enriched sequence motifs from data."""

    def __init__(self, reference_fasta, flank_size=2, min_sites=50):
        self.reference = reference_fasta
        self.flank_size = flank_size
        self.min_sites = min_sites

    def detect(self, bed_path, top_n=5, meth_thresh=0.7, unmeth_thresh=0.3):
        """Find motifs enriched in methylated vs unmethylated sites.

        Returns list of (motif, enrichment_ratio, p_meth, p_unmeth) sorted by enrichment.
        """
        from collections import Counter

        bed_path = Path(bed_path)
        ref_path = Path(self.reference)
        if not ref_path.exists():
            print("  [MotifDetector] Reference FASTA not found, skipping auto-detect")
            return []

        # Load reference into memory
        ref_seqs = {}
        chrom = None
        lines = []
        with open(ref_path) as f:
            for line in f:
                if line.startswith(">"):
                    if chrom and lines:
                        ref_seqs[chrom] = "".join(lines).upper()
                    chrom = line[1:].strip().split()[0]
                    lines = []
                else:
                    lines.append(line.strip())
            if chrom and lines:
                ref_seqs[chrom] = "".join(lines).upper()

        k = 1 + 2 * self.flank_size  # total motif length
        meth_counter = Counter()
        unmeth_counter = Counter()

        meth_count = 0
        unmeth_count = 0

        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                chrom, start_str = parts[0], parts[1]
                try:
                    pos = int(start_str)
                except ValueError:
                    continue
                meth, unmeth = int(parts[4]), int(parts[5])
                total = meth + unmeth
                if total == 0:
                    continue
                level = meth / total

                ref_seq = ref_seqs.get(chrom)
                if ref_seq is None:
                    continue

                start = pos - self.flank_size
                end = pos + self.flank_size + 1
                if start < 0 or end > len(ref_seq):
                    continue

                context = ref_seq[start:end]
                if len(context) != k:
                    continue

                if level >= meth_thresh:
                    meth_counter[context] += 1
                    meth_count += 1
                elif level <= unmeth_thresh:
                    unmeth_counter[context] += 1
                    unmeth_count += 1

        if meth_count < self.min_sites or unmeth_count < self.min_sites:
            print(f"  [MotifDetector] Too few sites (meth={meth_count}, unmeth={unmeth_count}), "
                  f"skipping detection")
            return []

        # Compute enrichment
        candidates = set(meth_counter.keys()) | set(unmeth_counter.keys())
        scored = []
        for ctx in candidates:
            m = meth_counter.get(ctx, 0)
            u = unmeth_counter.get(ctx, 0)
            p_meth = m / meth_count
            p_unmeth = u / unmeth_count if unmeth_count > 0 else 0
            ratio = (p_meth / p_unmeth) if p_unmeth > 0 else float("inf")
            if p_meth >= 0.01:   # at least 1% occurrence
                scored.append((ctx, ratio, p_meth, p_unmeth))

        scored.sort(key=lambda x: -x[1])
        motifs = [(ctx, f"{r:.2f}", f"{p:.3f}", f"{q:.3f}")
                  for ctx, r, p, q in scored[:top_n]]
        return motifs

    @staticmethod
    def format_motifs(motifs):
        """Pretty-print detected motifs."""
        lines = ["  [MotifDetector] Top motifs (methylated enrichment):"]
        for ctx, ratio, p_meth, p_unmeth in motifs:
            lines.append(f"    {ctx}  enrichment={ratio}x  "
                         f"meth_frac={p_meth}  unmeth_frac={p_unmeth}")
        return "\n".join(lines)


class ContextCorrector:
    """Motif-aware Bayesian methylation level correction.

    For each methylation site, extract the flanking sequence motif,
    look up the motif-specific prior, and compute the posterior
    methylation level using a Beta-Binomial model.
    """

    def __init__(self, reference_fasta, motifs=None,
                 prior_strength=10, flank_size=2):
        self.reference = reference_fasta
        self.motifs = motifs or []
        self.prior_strength = prior_strength
        self.flank_size = flank_size
        self._motif_priors = {}

    def fit_priors(self, bed_path):
        """Estimate motif-specific prior methylation rates from data."""
        from collections import Counter

        bed_path = Path(bed_path)
        ref_path = Path(self.reference)
        if not ref_path.exists():
            print("  [ContextCorrector] Reference not found, falling back to flat prior")
            self._motif_priors = {}
            return

        # Load reference
        ref_seqs = {}
        chrom = None
        lines = []
        with open(ref_path) as f:
            for line in f:
                if line.startswith(">"):
                    if chrom and lines:
                        ref_seqs[chrom] = "".join(lines).upper()
                    chrom = line[1:].strip().split()[0]
                    lines = []
                else:
                    lines.append(line.strip())
            if chrom and lines:
                ref_seqs[chrom] = "".join(lines).upper()

        # Compile motif patterns
        patterns = []
        for m in self.motifs:
            try:
                patterns.append((m, re.compile(_iupac_to_regex(m))))
            except re.error:
                print(f"  [ContextCorrector] Skipping invalid motif: {m}")

        if not patterns:
            self._motif_priors = {}
            return

        # Classify each site by motif context
        k = 1 + 2 * self.flank_size
        motif_data = defaultdict(lambda: {"meth": 0, "unmeth": 0, "count": 0})

        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                chrom, start_str = parts[0], parts[1]
                try:
                    pos = int(start_str)
                except ValueError:
                    continue
                meth, unmeth = int(parts[4]), int(parts[5])
                total = meth + unmeth
                if total == 0:
                    continue

                ref_seq = ref_seqs.get(chrom)
                if ref_seq is None:
                    continue

                start = pos - self.flank_size
                end = pos + self.flank_size + 1
                if start < 0 or end > len(ref_seq):
                    continue

                context = ref_seq[start:end].upper()
                if len(context) != k:
                    continue

                matched = False
                for motif_name, regex in patterns:
                    if regex.search(context):
                        d = motif_data[motif_name]
                        d["meth"] += meth
                        d["unmeth"] += unmeth
                        d["count"] += 1
                        matched = True
                        break
                if not matched:
                    d = motif_data["_other"]
                    d["meth"] += meth
                    d["unmeth"] += unmeth
                    d["count"] += 1

        # Compute prior rates
        for motif_name, d in motif_data.items():
            total = d["meth"] + d["unmeth"]
            if total > 0:
                rate = d["meth"] / total
                self._motif_priors[motif_name] = (rate, self.prior_strength)

        # Print summary
        print(f"  [ContextCorrector] Motif priors estimated from {bed_path.name}:")
        for name, (rate, strength) in sorted(self._motif_priors.items()):
            print(f"    {name}: prior_rate={rate:.4f} (pseudo-count={strength})")

    def correct_bed(self, bed_path, output_path, min_depth=5):
        """Apply context-aware Bayesian correction to methylation BED."""
        bed_path = Path(bed_path)
        output_path = Path(output_path)

        ref_path = Path(self.reference)
        if not ref_path.exists():
            print("  [ContextCorrector] Reference not found, falling back to ratio")
            return self._fallback_correct(bed_path, output_path, min_depth)

        # Load reference
        ref_seqs = {}
        chrom = None
        lines = []
        with open(ref_path) as f:
            for line in f:
                if line.startswith(">"):
                    if chrom and lines:
                        ref_seqs[chrom] = "".join(lines).upper()
                    chrom = line[1:].strip().split()[0]
                    lines = []
                else:
                    lines.append(line.strip())
            if chrom and lines:
                ref_seqs[chrom] = "".join(lines).upper()

        # Compile motif regexes
        patterns = []
        for m in self.motifs:
            try:
                patterns.append((m, re.compile(_iupac_to_regex(m))))
            except re.error:
                pass

        k = 1 + 2 * self.flank_size
        corrected_cnt = 0
        total = 0

        with open(bed_path) as fin, open(output_path, "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    fout.write(line)
                    continue

                chrom, start_str = parts[0], parts[1]
                try:
                    pos = int(start_str)
                except ValueError:
                    fout.write(line)
                    continue

                meth = int(parts[4])
                unmeth = int(parts[5])
                depth = meth + unmeth
                if depth < min_depth:
                    fout.write(line)
                    continue

                total += 1

                # Determine motif context
                ref_seq = ref_seqs.get(chrom)
                if ref_seq is None:
                    raw_level = meth / depth
                    fout.write(f"{chrom}\t{pos}\t{pos+1}\t{raw_level:.4f}\t"
                               f"{meth}\t{unmeth}\t{depth}\n")
                    continue

                start = pos - self.flank_size
                end = pos + self.flank_size + 1
                if start < 0 or end > len(ref_seq):
                    raw_level = meth / depth
                    fout.write(f"{chrom}\t{pos}\t{pos+1}\t{raw_level:.4f}\t"
                               f"{meth}\t{unmeth}\t{depth}\n")
                    continue

                context = ref_seq[start:end].upper()
                if len(context) != k:
                    raw_level = meth / depth
                    fout.write(f"{chrom}\t{pos}\t{pos+1}\t{raw_level:.4f}\t"
                               f"{meth}\t{unmeth}\t{depth}\n")
                    continue

                # Find matching motif
                matched_motif = "_other"
                for motif_name, regex in patterns:
                    if regex.search(context):
                        matched_motif = motif_name
                        break

                # Apply Beta-Binomial posterior
                prior_rate, prior_strength = self._motif_priors.get(
                    matched_motif, (0.5, 2)
                )

                alpha = prior_rate * prior_strength + meth
                beta = (1 - prior_rate) * prior_strength + unmeth
                posterior_mean = alpha / (alpha + beta) if (alpha + beta) > 0 else 0

                if abs(posterior_mean - meth / depth) > 0.01:
                    corrected_cnt += 1

                fout.write(f"{chrom}\t{pos}\t{pos+1}\t{posterior_mean:.4f}\t"
                           f"{meth}\t{unmeth}\t{depth}\n")

        print(f"  [ContextCorrector] Corrected {corrected_cnt}/{total} sites "
               f"(min_depth={min_depth})")
        return str(output_path)

    def _fallback_correct(self, bed_path, output_path, min_depth):
        """Fallback: simple ratio-based (same as consensus)."""
        bed_path = Path(bed_path)
        output_path = Path(output_path)
        with open(bed_path) as fin, open(output_path, "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                if len(parts) >= 7:
                    meth, unmeth = int(parts[4]), int(parts[5])
                    depth = meth + unmeth
                    if depth >= min_depth:
                        level = meth / depth
                        fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t"
                                   f"{level:.4f}\t{meth}\t{unmeth}\t{depth}\n")
                    else:
                        fout.write(line)
                else:
                    fout.write(line)
        return str(output_path)


class ErrorCorrector:
    """Apply error detection and correction to aligned reads."""

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.detector = None  # lazy-init ReadErrorDetector for "model" strategy

    def correct_bed(self, bed_path, conversion, reference=None, min_depth=5,
                    min_agree_ratio=0.7):
        """Dispatch to appropriate BED correction strategy.

        Uses ContextCorrector if 'context' strategy is enabled and
        reference is available; otherwise falls back to consensus.
        """
        strategies = self.config.get("correction", "strategies", default=["consensus"])

        if "context" in strategies:
            read_length = self.config.get("correction", "read_length", default=150)
            algo, rec_flank, rec_prior = AlgorithmSelector.select(read_length)
            if algo == "simple":
                print(f"  [ErrorCorrector] read_length={read_length}bp â†?simple correction")
                return self.consensus_correct_bed(bed_path, min_depth=min_depth,
                                                  min_agree_ratio=min_agree_ratio)

            motifs = self.config.get("correction", "motifs", default=[])
            auto_detect = self.config.get("correction", "auto_detect_motifs", default=True)
            flank = self.config.get("correction", "motif_flank", default=rec_flank)
            prior_n = self.config.get("correction", "prior_strength", default=rec_prior)

            # If user didn't explicitly set flank/prior, use the recommended values
            if self.config.get("correction", "motif_flank") is None:
                flank = rec_flank
            if self.config.get("correction", "prior_strength") is None:
                prior_n = rec_prior

            ref = reference or self.config.get("correction", "reference")
            if ref and Path(ref).exists():
                corrector = ContextCorrector(
                    reference_fasta=ref,
                    motifs=list(motifs) if motifs else [],
                    prior_strength=prior_n,
                    flank_size=flank,
                )

                if not motifs and auto_detect:
                    detector = MotifDetector(ref, flank_size=flank)
                    detected = detector.detect(bed_path)
                    if detected:
                        print(detector.format_motifs(detected))
                        corrector.motifs = [m[0] for m in detected]
                    else:
                        print("  [ErrorCorrector] No motifs auto-detected, "
                              "using flat prior")

                corrector.fit_priors(bed_path)
                base_name = Path(bed_path).stem
                out_path = self.output_dir / f"{base_name}_corrected.bed"
                return corrector.correct_bed(bed_path, out_path, min_depth=min_depth)
            else:
                print("  [ErrorCorrector] Reference required for context correction, "
                      "falling back to consensus")

        return self.consensus_correct_bed(bed_path, min_depth=min_depth,
                                          min_agree_ratio=min_agree_ratio)

    def correct_sam(self, sam_path, conversion, strategies=None):
        """Read SAM, apply correction strategies, write corrected SAM."""
        if strategies is None:
            strategies = ["mq", "clip"]

        # Model-based correction handled as a separate pass
        if "model" in strategies:
            if self.detector is not None:
                sam_path, _ = self.detector.filter_sam(
                    str(sam_path), conversion,
                    threshold=self.config.get("correction", "model_threshold", default=0.5),
                    output_sam=self.output_dir / f"{Path(sam_path).stem.replace('.sam', '')}_model_filtered.sam"
                )
            else:
                print("  [ErrorCorrector] Model not trained, skipping model strategy")

        sam_path = Path(sam_path)
        base_name = sam_path.stem.replace(".sam", "")
        corrected_path = self.output_dir / f"{base_name}_corrected.sam"
        rejected_path = self.output_dir / f"{base_name}_rejected.txt"
        unconverted_fastq_path = self.output_dir / f"{base_name}_unconverted.fastq"

        stats = {"total": 0, "kept": 0, "rejected": {"mq": 0, "clip": 0, "pair": 0, "unconverted": 0}}
        rejected_records = []

        mq_threshold = self.config.get("correction", "min_mq", default=20)
        max_clip_pct = self.config.get("correction", "max_clip_pct", default=50)
        unconverted_max = self.config.get("correction", "max_unconverted", default=3)
        target_base = conversion.target_base
        converted_base = conversion.converted_base

        with open(sam_path) as fin, open(corrected_path, "w") as fout:
            if "unconverted" in strategies:
                fout_unconv = open(unconverted_fastq_path, "w")
            else:
                fout_unconv = None

            for line in fin:
                if line.startswith("@"):
                    fout.write(line)
                    continue

                stats["total"] += 1
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    fout.write(line)
                    stats["kept"] += 1
                    continue

                flag = int(parts[1])
                if flag & 0x4:
                    fout.write(line)
                    stats["kept"] += 1
                    continue

                reject_reason = None

                if "mq" in strategies:
                    try:
                        mq = int(parts[4])
                        if mq < mq_threshold:
                            reject_reason = "mq"
                            stats["rejected"]["mq"] += 1
                    except (ValueError, IndexError):
                        pass

                if reject_reason is None and "clip" in strategies:
                    cigar = parts[5]
                    clip_bases = self._count_clipped_bases(cigar)
                    seq_len = len(parts[9])
                    if seq_len > 0 and (clip_bases / seq_len * 100) > max_clip_pct:
                        reject_reason = "clip"
                        stats["rejected"]["clip"] += 1

                if reject_reason is None and "pair" in strategies:
                    if self._is_discordant_pair(parts):
                        reject_reason = "pair"
                        stats["rejected"]["pair"] += 1

                if reject_reason is None and "unconverted" in strategies:
                    seq = parts[9].upper()
                    target_count = seq.count(target_base)
                    if target_count > unconverted_max:
                        reject_reason = "unconverted"
                        stats["rejected"]["unconverted"] += 1
                        if fout_unconv:
                            fout_unconv.write(f"@{parts[0]}\n{seq}\n+\n{parts[10]}\n")

                if reject_reason:
                    rejected_records.append((parts[0], reject_reason))
                else:
                    fout.write(line)
                    stats["kept"] += 1

            if fout_unconv:
                fout_unconv.close()

        with open(rejected_path, "w") as f:
            f.write("read_name\treason\n")
            for name, reason in rejected_records:
                f.write(f"{name}\t{reason}\n")

        print(f"  [ErrorCorrector] Total: {stats['total']}, Kept: {stats['kept']}, "
              f"Rejected: {sum(stats['rejected'].values())}")
        for reason, count in stats["rejected"].items():
            if count:
                print(f"    - {reason}: {count}")

        return str(corrected_path)

    def train_model(self, sam_path, ground_truth_bed, conversion):
        """Train the read error detection model."""
        from .model_correction import ReadErrorDetector
        print("  [ErrorCorrector] Training read error detection model...")
        detector = ReadErrorDetector()
        detector.train(sam_path, ground_truth_bed, conversion)
        self.detector = detector
        # Save model for reuse
        model_path = self.output_dir / "read_error_model.json"
        detector.save(model_path)
        print(f"  [ErrorCorrector] Model saved to {model_path}")
        return self.detector

    def consensus_correct_bed(self, bed_path, min_depth=5, min_agree_ratio=0.7):
        """Apply simple consensus correction (hard threshold aggregation)."""
        bed_path = Path(bed_path)
        base_name = bed_path.stem
        out_path = self.output_dir / f"{base_name}_consensus.bed"

        sites = defaultdict(list)
        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    key = (parts[0], int(parts[1]))
                    meth = int(parts[4])
                    unmeth = int(parts[5])
                    total = meth + unmeth
                    sites[key].append((meth, unmeth, total))

        corrected = 0
        with open(out_path, "w") as f:
            for (chrom, pos), calls in sorted(sites.items()):
                total_meth = sum(c[0] for c in calls)
                total_unmeth = sum(c[1] for c in calls)
                total_depth = total_meth + total_unmeth

                if total_depth < min_depth:
                    continue

                meth_ratio = total_meth / total_depth if total_depth > 0 else 0
                level = meth_ratio

                if meth_ratio >= min_agree_ratio:
                    level = 1.0
                    corrected += 1
                elif meth_ratio <= (1 - min_agree_ratio):
                    level = 0.0
                    corrected += 1

                f.write(f"{chrom}\t{pos}\t{pos+1}\t{level:.4f}\t"
                        f"{total_meth}\t{total_unmeth}\t{total_depth}\n")

        print(f"  [ErrorCorrector] Consensus correction: {corrected} sites adjusted "
              f"(depth >= {min_depth}, agree >= {min_agree_ratio})")
        return str(out_path)

    @staticmethod
    def _count_clipped_bases(cigar):
        clip = 0
        for m in re.finditer(r'(\d+)([SH])', cigar):
            clip += int(m.group(1))
        return clip

    @staticmethod
    def _is_discordant_pair(parts):
        flag = int(parts[1])
        if not (flag & 0x1):
            return False
        if flag & 0x4 and flag & 0x8:
            return False
        tlen = int(parts[8]) if len(parts) > 8 else 0
        rnext = parts[6]
        if rnext not in ("*", "="):
            return True
        if abs(tlen) > 1000:
            return True
        return False
