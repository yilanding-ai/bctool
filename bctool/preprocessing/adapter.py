"""
Adapter trimming module for BCTool.

Supports:
- Built-in adapter trimming (no external dependencies)
- Cutadapt integration (when available)
- Auto-detection of adapter sequences
- Single-end and paired-end
- C-to-T / A-to-G bisulfite-aware trimming
"""

import re
import gzip
import shutil
import subprocess
from pathlib import Path
from collections import Counter
from ..config import Config


# Common adapter sequences
ADAPTER_DB = {
    "illumina_universal": "AGATCGGAAGAGC",
    "illumina_small_rna": "TGGAATTCTCGG",
    "nextera": "CTGTCTCTTATA",
    "polya": "AAAAAAAAAAAA",
    "polyt": "TTTTTTTTTTTT",
    "truseq_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "truseq_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "nebnext": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "solexa": "TCGTATGCCGTCTTCTGCTTG",
}


class AdapterTrimmer:
    """Built-in adapter trimmer with auto-detection and cutadapt support."""

    def __init__(self, config: Config, output_dir: str):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.conversion = config.conversion
        self.cutadapt_available = shutil.which("cutadapt") is not None

    def auto_detect_adapters(self, read1_fastq, num_reads=10000):
        """Auto-detect adapter sequences by finding overrepresented kmers in read ends."""
        seqs_r1_start = []
        seqs_r1_end = []
        seqs_r2_start = []
        seqs_r2_end = []

        def _sample(path, start_list, end_list):
            if not path or not Path(path).exists():
                return
            opener = gzip.open if str(path).endswith(".gz") else open
            mode = "rt" if str(path).endswith(".gz") else "r"
            with opener(path, mode) as f:
                for i, line in enumerate(f):
                    if i >= num_reads * 4:
                        break
                    if i % 4 == 1:
                        seq = line.strip()
                        if len(seq) >= 10:
                            start_list.append(seq[:8])
                            end_list.append(seq[-12:])

        _sample(read1_fastq, seqs_r1_start, seqs_r1_end)

        read2 = self.config.read2
        if read2:
            _sample(read2, seqs_r2_start, seqs_r2_end)

        adapters = {}

        # Detect adapter at 3' end (most common location)
        for label, end_list, start_list in [
            ("r1_3prime", seqs_r1_end, seqs_r1_start),
        ]:
            counter = Counter(end_list)
            for seq, count in counter.most_common(20):
                if count < len(end_list) * 0.01:
                    continue
                for ad_name, ad_seq in ADAPTER_DB.items():
                    ad_len = min(len(ad_seq), 12)
                    # Check if adapter suffix matches read end
                    for i in range(len(ad_seq) - 5):
                        fragment = ad_seq[i:i + 8]
                        if fragment in seq or seq in fragment:
                            adapters[label] = ad_seq
                            print(f"  [Adapter] Detected {ad_name} at {label}: {ad_seq}")
                            break
                    if label in adapters:
                        break

        return adapters

    def trim(self, read1, read2=None, adapters=None):
        """Run adapter trimming. Uses cutadapt if available, else built-in."""
        if adapters is None:
            adapters = self._get_default_adapters()

        if self.cutadapt_available:
            return self._trim_cutadapt(read1, read2, adapters)
        else:
            return self._trim_builtin(read1, read2, adapters)

    def _get_default_adapters(self):
        """Get default adapters from config or use standard Illumina."""
        adapters = {}

        adapters["r1_3prime"] = self.config.get("adapter", "r1_3prime",
                                                  default=ADAPTER_DB["illumina_universal"])
        if self.config.mode == "pe":
            adapters["r2_3prime"] = self.config.get("adapter", "r2_3prime",
                                                     default=ADAPTER_DB["illumina_universal"])
        return adapters

    def _trim_cutadapt(self, read1, read2, adapters):
        """Use cutadapt for adapter trimming."""
        out_r1 = str(self.output_dir / "trimmed_R1.fastq")
        out_r2 = str(self.output_dir / "trimmed_R2.fastq") if read2 else None
        report = str(self.output_dir / "cutadapt_report.txt")

        cmd = ["cutadapt", "-j", str(self.config.threads)]

        for key, seq in adapters.items():
            if key.endswith("_3prime") or key == "r1_3prime":
                cmd.extend(["-a", seq])
            elif key.endswith("_5prime"):
                cmd.extend(["-g", seq])

        cmd.extend(["-o", out_r1])
        if out_r2:
            cmd.extend(["-p", out_r2])
        cmd.extend([read1])
        if read2:
            cmd.append(read2)

        cmd_str = " ".join(cmd) + f" > {report}"
        print(f"  [cutadapt] Running: cutadapt (report: {report})")
        result = subprocess.run(cmd, capture_output=True, text=True)

        with open(report, "w") as f:
            f.write(result.stdout)
            f.write("\n--- stderr ---\n")
            f.write(result.stderr)

        return out_r1, out_r2

    def _trim_builtin(self, read1, read2, adapters):
        """Built-in adapter trimming (no external dependencies)."""
        out_r1 = str(self.output_dir / "trimmed_R1.fastq")
        out_r2 = str(self.output_dir / "trimmed_R2.fastq") if read2 else None

        r1_adapter = adapters.get("r1_3prime", ADAPTER_DB["illumina_universal"])
        r2_adapter = adapters.get("r2_3prime", r1_adapter)

        stats = {"r1_trimmed": 0, "r2_trimmed": 0, "r1_kept": 0, "r2_kept": 0}

        def _trim_reads(in_path, out_path, adapter):
            if not in_path or not Path(in_path).exists():
                return 0
            trimmed = 0
            kept = 0
            opener = gzip.open if str(in_path).endswith(".gz") else open
            mode = "rt" if str(in_path).endswith(".gz") else "r"
            with opener(in_path, mode) as fin, open(out_path, "w") as fout:
                while True:
                    header = fin.readline()
                    if not header:
                        break
                    seq = fin.readline().rstrip()
                    plus = fin.readline()
                    qual = fin.readline().rstrip()

                    trimmed_seq, trimmed_qual = self._trim_adapter(seq, qual, adapter)
                    if trimmed_seq != seq:
                        trimmed += 1
                    kept += 1

                    fout.write(f"{header}{trimmed_seq}\n{plus}{trimmed_qual}\n")
            return trimmed

        stats["r1_trimmed"] = _trim_reads(read1, out_r1, r1_adapter)
        stats["r1_kept"] = self._count_reads(out_r1)

        if read2:
            stats["r2_trimmed"] = _trim_reads(read2, out_r2, r2_adapter)
            stats["r2_kept"] = self._count_reads(out_r2)

        total = stats["r1_trimmed"] + stats.get("r2_trimmed", 0)
        print(f"  [Trimmer] Built-in: trimmed {total} adapters "
              f"(R1: {stats['r1_trimmed']}, R2: {stats.get('r2_trimmed', 0)})")
        return out_r1, out_r2

    def _trim_adapter(self, seq, qual, adapter, min_overlap=6, max_error=0.1):
        """Trim adapter from 3' end of sequence using semi-global alignment."""
        if len(seq) < min_overlap:
            return seq, qual

        # For bisulfite-converted reads, also search with conversion-aware adapter
        adapters_to_check = [adapter]
        from ..methylation.converter import ConversionType
        conv = ConversionType(self.conversion)
        if conv.target_base in "ACGT" and conv.converted_base in "ACGT":
            adapters_to_check.append(
                adapter.replace(conv.target_base, conv.converted_base)
            )

        best_pos = len(seq)
        for ad in adapters_to_check:
            for i in range(len(seq) - min_overlap + 1):
                fragment = seq[i:]
                overlap = min(len(fragment), len(ad))
                if overlap < min_overlap:
                    continue
                matches = sum(1 for a, b in zip(fragment[:overlap], ad[:overlap])
                              if a == b)
                error_rate = 1 - matches / overlap
                if error_rate <= max_error and matches >= min_overlap:
                    if i < best_pos:
                        best_pos = i
                    break

        if best_pos < len(seq):
            return seq[:best_pos], qual[:best_pos]
        return seq, qual

    @staticmethod
    def _count_reads(fastq_path):
        if not Path(fastq_path).exists():
            return 0
        n = 0
        with open(fastq_path) as f:
            for _ in f:
                n += 1
        return n // 4

    @staticmethod
    def list_known_adapters():
        """Print known adapter sequences."""
        print("Known adapter sequences:")
        print(f"  {'Name':<25} {'Sequence'}")
        print("-" * 50)
        for name, seq in ADAPTER_DB.items():
            print(f"  {name:<25} {seq}")
