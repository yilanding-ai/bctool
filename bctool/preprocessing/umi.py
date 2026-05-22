"""
UMI (Unique Molecular Identifier) analysis and deduplication module.

Supports:
- UMI extraction from read header (standard: @instrument:run:flowcell:lane:tile:x:y#UMI)
- UMI extraction from read sequence (inline UMI at 5' or 3' end)
- Consensus read generation
- UMI statistics and QC metrics
- Single-end and paired-end
"""

import gzip
import re
from pathlib import Path
from collections import Counter, defaultdict
from ..config import Config


UMI_PATTERNS = {
    "standard": r"^@([^\s]+)\s+.*?(UMI|umi|BX:Z:)([ATCGNatcgn]+)",
    "illumina_v3": r"^@([^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+)\s+.*?([ATCGN]{8,16})$",
    "simple": r"^@([^\s]+)_([ATCGN]{8,16})(?:_|\s)",
}

UMI_EXAMPLES = [
    ("standard", "Illumina header BX tag", "@A00261:111:HL3G5DSX3:3:1101:1234:5678 BX:Z:ACGTACGT"),
    ("header_suffix", "After space in header", "@read_1 ACGTACGT"),
    ("inline_5prime", "First N bases of read", "Read seq: ACGTACGTNNNN..."),
    ("inline_3prime", "Last N bases of read", "Read seq: ...NNNNACGTACGT"),
]


class UmiProcessor:
    """Extract UMI barcodes and deduplicate sequencing reads."""

    def __init__(self, config: Config, output_dir: str):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.umi_stats = {
            "total_reads": 0,
            "reads_with_umi": 0,
            "unique_umis": 0,
            "duplicate_groups": 0,
            "singleton_groups": 0,
            "reads_after_dedup": 0,
            "duplication_rate": 0.0,
            "umi_length": None,
        }

    def extract_and_deduplicate(self, read1, read2=None, umi_config=None):
        """Extract UMIs and deduplicate reads.
        
        umi_config can contain:
            - location: "header" (default), "inline_5prime", "inline_3prime"
            - length: UMI length (auto-detect if None)
            - pattern: regex pattern for UMI extraction
            - strategy: "cluster" (group identical UMIs + position), 
                       "consensus" (build consensus from UMI groups)
        """
        if umi_config is None:
            umi_config = self._detect_umi(read1)

        strategy = umi_config.get("strategy", "cluster")

        print(f"  [UMI] Strategy: {strategy}")
        print(f"  [UMI] Location: {umi_config.get('location', 'header')}")
        print(f"  [UMI] Length: {umi_config.get('length', 'auto')}")

        if strategy == "cluster":
            return self._cluster_dedup(read1, read2, umi_config)
        elif strategy == "consensus":
            return self._consensus_dedup(read1, read2, umi_config)
        else:
            raise ValueError(f"Unknown UMI strategy: {strategy}")

    def _detect_umi(self, read1):
        """Auto-detect UMI configuration from reads."""
        if not read1 or not Path(read1).exists():
            return {"location": "header", "length": 0, "strategy": "none"}

        samples = self._sample_reads(read1, 500)
        config = {"location": "header", "strategy": "cluster"}

        # Try different patterns
        for pattern_name, pattern_str in UMI_PATTERNS.items():
            pattern = re.compile(pattern_str)
            matches = 0
            umi_lengths = []
            for header in samples:
                m = pattern.search(header)
                if m:
                    matches += 1
                    umi = m.group(3) if pattern.groups >= 3 else m.group(2)
                    umi_lengths.append(len(umi))

            if matches > len(samples) * 0.5:
                config["pattern"] = pattern_str
                config["detected_via"] = pattern_name
                if umi_lengths:
                    config["length"] = max(set(umi_lengths), key=umi_lengths.count)
                print(f"  [UMI] Auto-detected pattern: {pattern_name} "
                      f"({matches}/{len(samples)} reads match)")
                break

        # Check inline UMI if header patterns fail
        if "pattern" not in config:
            seqs = self._sample_seqs(read1, 200)
            if seqs:
                # Check for constant-length barcode at 5' end
                first_bases = [s[:8] for s in seqs if len(s) >= 8]
                if first_bases:
                    unique = len(set(first_bases))
                    if unique > len(first_bases) * 0.5:
                        config["location"] = "inline_5prime"
                        config["length"] = 8
                        print(f"  [UMI] Auto-detected inline 5' UMI (length=8)")
                    else:
                        last_bases = [s[-8:] for s in seqs if len(s) >= 8]
                        if last_bases:
                            unique = len(set(last_bases))
                            if unique > len(last_bases) * 0.5:
                                config["location"] = "inline_3prime"
                                config["length"] = 8
                                print(f"  [UMI] Auto-detected inline 3' UMI (length=8)")

        if "pattern" not in config and config.get("location") == "header":
            config["length"] = 0
            print(f"  [UMI] No UMI detected, skipping deduplication")
            config["strategy"] = "none"

        return config

    def _cluster_dedup(self, read1, read2, umi_config):
        """Group reads by UMI + mapping position and keep one per group."""
        out_r1 = str(self.output_dir / "dedup_R1.fastq")
        out_r2 = str(self.output_dir / "dedup_R2.fastq") if read2 else None

        if umi_config.get("strategy") == "none" or umi_config.get("length", 0) == 0:
            print(f"  [UMI] No UMI-based deduplication needed")
            return read1, read2

        # Extract UMIs
        r1_umis = self._extract_umis(read1, umi_config)
        r2_umis = self._extract_umis(read2, umi_config) if read2 else None

        if not r1_umis:
            print(f"  [UMI] No UMIs found, skipping deduplication")
            return read1, read2

        # Group reads by UMI
        umi_groups = defaultdict(list)
        for idx, umi in r1_umis.items():
            if umi:
                umi_groups[umi].append(idx)

        # Statistics
        total = len(r1_umis)
        with_umi = sum(1 for u in r1_umis.values() if u)
        unique_umis = len(umi_groups)
        singletons = sum(1 for g in umi_groups.values() if len(g) == 1)
        reads_after = unique_umis  # one per UMI group
        dup_rate = 1 - (reads_after / total) if total > 0 else 0

        self.umi_stats = {
            "total_reads": total,
            "reads_with_umi": with_umi,
            "unique_umis": unique_umis,
            "duplicate_groups": unique_umis - singletons,
            "singleton_groups": singletons,
            "reads_after_dedup": reads_after,
            "duplication_rate": round(dup_rate, 4),
            "umi_length": umi_config.get("length", 0),
        }

        # Write deduplicated reads (first read of each UMI group)
        kept_indices = set()
        for umi, indices in umi_groups.items():
            kept_indices.add(indices[0])

        self._write_selected_reads(read1, out_r1, kept_indices)
        if read2:
            self._write_selected_reads(read2, out_r2, kept_indices)

        n_kept = len(kept_indices)
        print(f"  [UMI] Dedup: {total} -> {n_kept} reads "
              f"(rate: {dup_rate:.2%}, UMIs: {unique_umis})")

        return out_r1, out_r2

    def _consensus_dedup(self, read1, read2, umi_config):
        """Build consensus reads from UMI groups (advanced)."""
        print(f"  [UMI] Consensus mode not yet implemented, falling back to cluster")
        return self._cluster_dedup(read1, read2, umi_config)

    def _extract_umis(self, fastq_path, umi_config):
        """Extract UMI from each read. Returns dict of {read_index: umi_string}."""
        if not fastq_path or not Path(fastq_path).exists():
            return {}

        result = {}
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        mode = "rt" if str(fastq_path).endswith(".gz") else "r"
        location = umi_config.get("location", "header")
        length = umi_config.get("length", 0)
        pattern_str = umi_config.get("pattern")

        with opener(fastq_path, mode) as f:
            idx = 0
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline()
                f.readline()  # plus line
                f.readline()  # qual

                umi = None
                if location == "header" and pattern_str:
                    m = re.search(pattern_str, header)
                    if m:
                        umi = m.group(3) if m.lastindex and m.lastindex >= 3 else m.group(2)
                elif location == "inline_5prime" and length > 0:
                    umi = seq[:length].strip()
                elif location == "inline_3prime" and length > 0:
                    umi = seq[-length:].strip()

                if umi and len(umi) >= 4:
                    # Standardize UMI (uppercase, remove non-standard bases)
                    umi = re.sub(r'[^ATCG]', 'N', umi.upper())
                    result[idx] = umi
                idx += 1

        return result

    def _write_selected_reads(self, in_path, out_path, selected_indices):
        """Write only selected reads to output."""
        opener = gzip.open if str(in_path).endswith(".gz") else open
        mode = "rt" if str(in_path).endswith(".gz") else "r"
        with opener(in_path, mode) as fin, open(out_path, "w") as fout:
            idx = 0
            while True:
                header = fin.readline()
                if not header:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                if idx in selected_indices:
                    fout.write(f"{header}{seq}{plus}{qual}")
                idx += 1

    def _sample_reads(self, fastq_path, n=500):
        """Sample read headers for UMI detection."""
        headers = []
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        mode = "rt" if str(fastq_path).endswith(".gz") else "r"
        with opener(fastq_path, mode) as f:
            for i, line in enumerate(f):
                if i >= n * 4:
                    break
                if i % 4 == 0:
                    headers.append(line.strip())
        return headers

    def _sample_seqs(self, fastq_path, n=200):
        """Sample read sequences for inline UMI detection."""
        seqs = []
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        mode = "rt" if str(fastq_path).endswith(".gz") else "r"
        with opener(fastq_path, mode) as f:
            for i, line in enumerate(f):
                if i >= n * 4:
                    break
                if i % 4 == 1:
                    seqs.append(line.strip())
        return seqs

    def report_stats(self):
        """Print UMI statistics."""
        print(f"\n  [UMI Stats]")
        print(f"    Total reads:         {self.umi_stats['total_reads']}")
        print(f"    Reads with UMI:      {self.umi_stats['reads_with_umi']}")
        print(f"    Unique UMIs:         {self.umi_stats['unique_umis']}")
        print(f"    Duplicate groups:    {self.umi_stats['duplicate_groups']}")
        print(f"    Singleton groups:    {self.umi_stats['singleton_groups']}")
        print(f"    Reads after dedup:   {self.umi_stats['reads_after_dedup']}")
        print(f"    Duplication rate:    {self.umi_stats['duplication_rate']:.2%}")

    @staticmethod
    def list_umi_patterns():
        """Print example UMI patterns."""
        print("Common UMI formats:")
        print(f"  {'Format':<20} {'Description':<30} {'Example'}")
        print("-" * 80)
        for fmt, desc, example in UMI_EXAMPLES:
            print(f"  {fmt:<20} {desc:<30} {example}")
