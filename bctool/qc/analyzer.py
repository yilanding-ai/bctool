"""
Quality Control analyzer for FASTQ sequencing data.

Provides FastQC-equivalent metrics:
- Per-base quality scores
- Per-base sequence content (A/T/G/C/N)
- Per-base GC content
- Per-sequence quality scores
- Per-sequence GC content
- Sequence length distribution
- N content per base
- Overrepresented sequences
- Adapter content detection
"""

import gzip
import math
import re
from pathlib import Path
from collections import Counter, defaultdict


class QCAnalyzer:
    """Analyze FASTQ files and compute quality metrics."""

    def __init__(self, max_reads=100000):
        self.max_reads = max_reads
        self.metrics = {}
        self._reset()

    def _reset(self):
        self.metrics = {
            "total_reads": 0,
            "total_bases": 0,
            "read_lengths": Counter(),
            "base_counts": Counter(),        # A, C, G, T, N
            "per_base_qual": defaultdict(list),
            "per_base_content": defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}),
            "per_base_gc": defaultdict(list),
            "per_seq_qual": [],
            "per_seq_gc": [],
            "n_content": defaultdict(int),
            "sequences": [],
            "gc_bins": Counter(),
            "qual_bins": Counter(),
        }

    def analyze(self, fastq_path, read2_path=None):
        """Run full QC analysis on FASTQ file(s)."""
        self._reset()
        path = Path(fastq_path)

        if not path.exists():
            print(f"  [QC] File not found: {fastq_path}")
            return self.metrics

        print(f"  [QC] Analyzing {path.name}...")
        self._analyze_single(path)
        self._compute_summary_stats()

        if read2_path and Path(read2_path).exists():
            print(f"  [QC] Analyzing {Path(read2_path).name}...")
            self._analyze_single(read2_path)
            self._compute_summary_stats()

        self._find_overrepresented()
        return self.metrics

    def _analyze_single(self, fastq_path):
        """Analyze a single FASTQ file."""
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        mode = "rt" if str(fastq_path).endswith(".gz") else "r"

        with opener(fastq_path, mode) as f:
            read_idx = 0
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().rstrip()
                plus = f.readline()
                qual = f.readline().rstrip()

                self._process_read(seq, qual)
                read_idx += 1
                if read_idx >= self.max_reads:
                    break

    def _process_read(self, seq, qual):
        """Process a single read pair."""
        seq_len = len(seq)
        self.metrics["read_lengths"][seq_len] += 1
        self.metrics["total_reads"] += 1
        self.metrics["total_bases"] += seq_len

        # Per-base metrics
        for i, (base, q) in enumerate(zip(seq, qual)):
            base = base.upper()
            if base not in "ACGTN":
                base = "N"

            self.metrics["base_counts"][base] += 1
            self.metrics["per_base_content"][i][base] += 1

            if base == "N":
                self.metrics["n_content"][i] += 1

            q_score = ord(q) - 33
            if 0 <= q_score <= 41:
                self.metrics["per_base_qual"][i].append(q_score)

        # Per-sequence metrics
        q_scores = [min(41, max(0, ord(q) - 33)) for q in qual[:seq_len]]
        mean_qual = sum(q_scores) / len(q_scores) if q_scores else 0
        self.metrics["per_seq_qual"].append(mean_qual)
        self.metrics["qual_bins"][int(mean_qual)] += 1

        gc_count = seq.upper().count("G") + seq.upper().count("C")
        gc_pct = gc_count / seq_len * 100 if seq_len > 0 else 0
        self.metrics["per_seq_gc"].append(gc_pct)
        self.metrics["gc_bins"][int(gc_pct)] += 1

        # Store for overrepresented detection
        self.metrics["sequences"].append(seq[:50])

    def _compute_summary_stats(self):
        """Compute summary statistics from raw metrics."""
        m = self.metrics
        base_total = sum(m["base_counts"].values()) or 1

        # Base composition percentages
        m["base_pct"] = {
            b: m["base_counts"].get(b, 0) / base_total * 100
            for b in "ACGTN"
        }

        # Per-base quality summary (mean, median, quartiles)
        m["per_base_qual_summary"] = {}
        for pos, vals in m["per_base_qual"].items():
            if not vals:
                continue
            sorted_vals = sorted(vals)
            n = len(sorted_vals)
            m["per_base_qual_summary"][pos] = {
                "mean": sum(sorted_vals) / n,
                "med": sorted_vals[n // 2],
                "q1": sorted_vals[n // 4],
                "q3": sorted_vals[3 * n // 4],
                "lower": sorted_vals[max(0, int(n * 0.1))],
                "upper": sorted_vals[min(n - 1, int(n * 0.9))],
                "min": sorted_vals[0],
                "max": sorted_vals[-1],
            }

        # Per-base GC content
        for pos, counts in m["per_base_content"].items():
            total = sum(counts.values()) or 1
            m["per_base_gc"][pos] = (counts["G"] + counts["C"]) / total * 100

        # Sequence length summary
        lengths = list(m["read_lengths"].elements())
        m["length_summary"] = {
            "min": min(lengths) if lengths else 0,
            "max": max(lengths) if lengths else 0,
            "mean": sum(lengths) / len(lengths) if lengths else 0,
            "n50": self._n50(lengths) if lengths else 0,
        }

        # GC content summary
        gc_vals = m["per_seq_gc"]
        m["gc_summary"] = {
            "mean": sum(gc_vals) / len(gc_vals) if gc_vals else 0,
            "min": min(gc_vals) if gc_vals else 0,
            "max": max(gc_vals) if gc_vals else 0,
        }

        # Quality summary
        q_vals = m["per_seq_qual"]
        m["qual_summary"] = {
            "mean": sum(q_vals) / len(q_vals) if q_vals else 0,
            "min": min(q_vals) if q_vals else 0,
            "max": max(q_vals) if q_vals else 0,
        }

    def _find_overrepresented(self, top_n=20):
        """Find overrepresented sequences (potential adapter contamination)."""
        seqs = self.metrics["sequences"]
        if not seqs:
            return

        counter = Counter(seqs)
        total = len(seqs)
        overrep = []
        for seq, count in counter.most_common(top_n):
            pct = count / total * 100
            if pct >= 0.1:
                overrep.append({
                    "sequence": seq,
                    "count": count,
                    "pct": round(pct, 2),
                })

        self.metrics["overrepresented"] = overrep

        # Adapter detection from overrepresented
        self._detect_adapters(overrep)

    def _detect_adapters(self, overrepresented):
        """Match overrepresented sequences against known adapters."""
        known_adapters = {
            "Illumina Universal": "AGATCGGAAGAGC",
            "Illumina Small RNA": "TGGAATTCTCGG",
            "Nextera": "CTGTCTCTTATA",
            "Truseq R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "Truseq R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        }

        detected = []
        for item in overrepresented:
            seq = item["sequence"]
            for name, ad_seq in known_adapters.items():
                # Check for adapter match in first/last 30bp
                overlap = max(0, min(len(seq), len(ad_seq)) - 5)
                if overlap < 8:
                    continue
                for start in range(len(ad_seq) - overlap + 1):
                    frag = ad_seq[start:start + overlap]
                    if frag in seq or seq[:len(frag)] == frag or seq[-len(frag):] == frag:
                        detected.append({
                            "adapter": name,
                            "sequence": ad_seq,
                            "match": item["sequence"],
                            "pct": item["pct"],
                        })
                        break
                if detected and detected[-1]["adapter"] == name:
                    break

        self.metrics["detected_adapters"] = detected

    def summary_text(self):
        """Generate QC summary text."""
        m = self.metrics
        lines = []
        lines.append("  ┌─────────────────────────────────────────────�?)
        lines.append("  �?          QC Summary Report                 �?)
        lines.append("  ├─────────────────────────────────────────────�?)
        lines.append(f"  �?Total reads:      {m['total_reads']:<12}        �?)
        lines.append(f"  �?Total bases:      {m['total_bases']:<12}        �?)

        if "length_summary" in m:
            ls = m["length_summary"]
            lines.append(f"  �?Read length:      {ls['mean']:.0f} (N50:{ls['n50']})        �?)

        if "base_pct" in m:
            bp = m["base_pct"]
            gc = bp.get("G", 0) + bp.get("C", 0)
            lines.append(f"  �?GC content:       {gc:.1f}%                    �?)

        if "qual_summary" in m:
            qs = m["qual_summary"]
            lines.append(f"  �?Mean qual:        {qs['mean']:.1f}                    �?)

        if "overrepresented" in m and m["overrepresented"]:
            lines.append(f"  �?Overrep seqs:     {len(m['overrepresented'])}                     �?)

        if "detected_adapters" in m and m["detected_adapters"]:
            for ad in m["detected_adapters"]:
                lines.append(f"  �?Adapter: {ad['adapter']:<20}    �?)

        lines.append("  └─────────────────────────────────────────────�?)
        return "\n".join(lines)

    def get_grade(self):
        """Return overall QC grade."""
        issues = 0
        m = self.metrics

        qual = m.get("qual_summary", {}).get("mean", 40)
        if qual < 20:
            issues += 3
        elif qual < 28:
            issues += 1

        gc = sum(m.get("base_pct", {}).get(b, 0) for b in "GC")
        if gc < 20 or gc > 80:
            issues += 2
        elif gc < 30 or gc > 70:
            issues += 1

        n_pct = m.get("base_pct", {}).get("N", 0)
        if n_pct > 5:
            issues += 2
        elif n_pct > 1:
            issues += 1

        if m.get("detected_adapters"):
            issues += 2

        if m.get("overrepresented") and len(m["overrepresented"]) > 5:
            issues += 1

        if issues == 0:
            return "PASS (green)"
        elif issues <= 2:
            return "WARNING (yellow)"
        else:
            return "FAIL (red)"

    @staticmethod
    def _n50(lengths):
        """Compute N50 of read lengths."""
        sorted_l = sorted(lengths, reverse=True)
        total = sum(sorted_l)
        half = total / 2
        running = 0
        for l in sorted_l:
            running += l
            if running >= half:
                return l
        return 0
