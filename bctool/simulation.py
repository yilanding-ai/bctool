import random
import gzip
from pathlib import Path
import numpy as np
from .methylation.converter import ConversionType


class Simulator:
    """Generate simulated bisulfite/converted sequencing data with known ground truth."""

    def __init__(self, genome_length=1_000_000, num_reads=50_000, read_length=100,
                 fragment_mean=300, fragment_sd=50, methylation_rate=0.7,
                 error_rate=0.01, conversion_rate=0.99,
                 conversion="ct", mode="pe", seed=42):
        self.genome_length = genome_length
        self.num_reads = num_reads
        self.read_length = read_length
        self.fragment_mean = fragment_mean
        self.fragment_sd = fragment_sd
        self.methylation_rate = methylation_rate
        self.error_rate = error_rate
        self.conversion_rate = conversion_rate
        self.conv = ConversionType(conversion)
        self.mode = mode
        self.seed = seed
        self.rng = random.Random(seed)

    def _generate_genome(self):
        bases = ["A", "C", "G", "T"]
        genome = "".join(self.rng.choices(bases, k=self.genome_length))
        return genome

    def _select_conversion_sites(self, genome):
        target_base = self.conv.target_base
        sites = [i for i, b in enumerate(genome) if b == target_base]
        methylated = set()
        for s in sites:
            if self.rng.random() < self.methylation_rate:
                methylated.add(s)
        return sites, methylated

    def _convert_base(self, base, is_methylated):
        conv = self.conv
        if base == conv.target_base:
            if not is_methylated:
                return conv.converted_base if self.rng.random() < self.conversion_rate else conv.target_base
            else:
                return conv.target_base if self.rng.random() < self.conversion_rate else conv.converted_base
        if base == conv.complement_target:
            if not is_methylated:
                return conv.complement_converted if self.rng.random() < self.conversion_rate else conv.complement_target
            else:
                return conv.complement_target if self.rng.random() < self.conversion_rate else conv.complement_converted
        return base

    def _add_errors(self, seq):
        bases = ["A", "C", "G", "T"]
        seq = list(seq)
        for i in range(len(seq)):
            if self.rng.random() < self.error_rate:
                orig = seq[i]
                choices = [b for b in bases if b != orig]
                seq[i] = self.rng.choice(choices)
        return "".join(seq)

    def _generate_qual(self, length, mean_qual=30):
        qual_chars = []
        for _ in range(length):
            q = min(41, max(2, int(self.rng.gauss(mean_qual, 5))))
            qual_chars.append(chr(q + 33))
        return "".join(qual_chars)

    def run(self, output_dir):
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        print("Generating reference genome...")
        genome = self._generate_genome()
        sites, methylated = self._select_conversion_sites(genome)

        genome_path = out / "reference.fa"
        with open(genome_path, "w") as f:
            f.write(">simulated_reference\n")
            for i in range(0, len(genome), 80):
                f.write(genome[i:i+80] + "\n")

        ground_truth = out / "ground_truth.bed"
        with open(ground_truth, "w") as f:
            for s in sites:
                status = "1" if s in methylated else "0"
                f.write(f"chr_sim\t{s}\t{s+1}\t{status}\n")
        print(f"Ground truth: {len(sites)} target sites, "
              f"{len(methylated)} methylated")

        print("Generating reads...")
        if self.mode == "pe":
            f1_records, f2_records = [], []
        else:
            f1_records = []

        genome_len = len(genome)
        for i in range(self.num_reads):
            if i % 10000 == 0 and i > 0:
                print(f"  {i}/{self.num_reads} reads generated")

            frag_len = max(self.read_length * 2 + 20,
                           int(self.rng.gauss(self.fragment_mean, self.fragment_sd)))
            start = self.rng.randint(0, max(0, genome_len - frag_len - 1))
            fragment = genome[start:start + frag_len]

            # convert fragment
            conv_frag = list(fragment)
            for j in range(len(conv_frag)):
                pos = start + j
                is_meth = pos in methylated
                conv_frag[j] = self._convert_base(conv_frag[j], is_meth)
            conv_frag = "".join(conv_frag)

            if self.mode == "pe":
                r1_seq = conv_frag[:self.read_length]
                r2_seq = self._reverse_complement(
                    conv_frag[-self.read_length:])
                r1_seq = self._add_errors(r1_seq)
                r2_seq = self._add_errors(r2_seq)
                qual1 = self._generate_qual(self.read_length)
                qual2 = self._generate_qual(self.read_length)
                r1_header = f"@sim_read_{i}/1"
                r2_header = f"@sim_read_{i}/2"
                f1_records.append((r1_header, r1_seq, qual1))
                f2_records.append((r2_header, r2_seq, qual2))
            else:
                read_start = self.rng.randint(0, max(0, frag_len - self.read_length))
                read_seq = conv_frag[read_start:read_start + self.read_length]
                read_seq = self._add_errors(read_seq)
                qual = self._generate_qual(self.read_length)
                header = f"@sim_read_{i}"
                f1_records.append((header, read_seq, qual))

        def _write_fastq(path, records):
            with open(path, "w") as f:
                for h, s, q in records:
                    f.write(f"{h}\n{s}\n+\n{q}\n")

        _write_fastq(out / "simulated_R1.fastq", f1_records)
        if self.mode == "pe":
            _write_fastq(out / "simulated_R2.fastq", f2_records)

        # write truth read positions
        print(f"Simulation complete. Files in {out}")
        return out

    @staticmethod
    def _reverse_complement(seq):
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(comp.get(b, b) for b in reversed(seq))
