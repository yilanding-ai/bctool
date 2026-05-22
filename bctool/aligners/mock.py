import random
from pathlib import Path
from .base import AlignerBase


class MockAligner(AlignerBase):
    """Built-in mock aligner that simulates alignment with configurable accuracy.
    
    This allows the pipeline to run end-to-end without any external tools.
    Each mock aligner simulates a different accuracy level to demonstrate
    the benchmarking and ranking capability.
    """

    name = "mock"
    binary_name = "__mock__"
    requires_index = False

    _instances = 0

    def __init__(self, config, output_dir, accuracy=0.85, name_suffix=""):
        super().__init__(config, output_dir)
        MockAligner._instances += 1
        if name_suffix:
            self.name = f"mock_{name_suffix}"
        else:
            self.name = f"mock_{MockAligner._instances}"
        self.aligner_dir = self.output_dir / self.name
        self.aligner_dir.mkdir(parents=True, exist_ok=True)
        self.accuracy = accuracy
        self.rng = random.Random(hash(self.name))

    def is_available(self):
        return True

    def build_index(self, reference_fa):
        pass

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        """Simulate alignment by copying reads and introducing noise."""
        out_sam = str(self.aligner_dir / "aligned.sam")
        self._generate_simulated_alignment(read1, read2, reference_fa, out_sam)
        return out_sam

    def call_methylation(self, sam_bam_path, reference_fa=None):
        """Simulate methylation calling with controlled accuracy."""
        from ..methylation.converter import ConversionType
        out_path = str(self.aligner_dir / "methylation.bed")

        conv = ConversionType(self.config.conversion)
        target = conv.target_base
        meth_char = conv.converted_base

        ref_path = Path(reference_fa) if reference_fa else None
        if ref_path and ref_path.exists():
            genome = self._read_genome(ref_path)
        else:
            genome = "N" * 1000000

        with open(out_path, "w") as f:
            for pos, base in enumerate(genome):
                if base.upper() == target:
                    noise = self.rng.random()
                    if noise < self.accuracy:
                        is_meth = self.rng.random() < 0.7
                    else:
                        is_meth = self.rng.random() < 0.3

                    level = 0.85 if is_meth else 0.12
                    meth = int(level * 10)
                    unmeth = 10 - meth
                    chrom = "chr_sim"
                    f.write(f"{chrom}\t{pos}\t{pos+1}\t{level:.4f}\t{meth}\t{unmeth}\t10\n")

        print(f"  [{self.name}] Generated simulated methylation (accuracy={self.accuracy})")
        return out_path

    def _generate_simulated_alignment(self, read1, read2, reference_fa, out_sam):
        ref_path = Path(reference_fa) if reference_fa else None
        genome_len = 0
        if ref_path and ref_path.exists():
            genome = self._read_genome(ref_path)
            genome_len = len(genome)
        else:
            genome = "N" * 1000000
            genome_len = 1000000

        with open(out_sam, "w") as f:
            f.write("@HD\tVN:1.6\tSO:coordinate\n")
            f.write(f"@SQ\tSN:chr_sim\tLN:{genome_len}\n")

            fastq_path = Path(read1)
            if not fastq_path.exists():
                print(f"  [{self.name}] Warning: {read1} not found, creating empty SAM")
                return

            with open(fastq_path) as fq:
                lines = fq.readlines()
            for i in range(0, len(lines), 4):
                if i + 3 >= len(lines):
                    break
                header = lines[i].strip().lstrip("@")
                seq = lines[i+1].strip()
                qual = lines[i+3].strip()
                pos = self.rng.randint(0, max(0, genome_len - len(seq)))
                flag = 0 if read2 is None else 99
                mapq = 60
                cigar = f"{len(seq)}M"
                f.write(f"{header}\t{flag}\tchr_sim\t{pos+1}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")

    @staticmethod
    def _read_genome(fa_path):
        seq = []
        with open(fa_path) as f:
            for line in f:
                if not line.startswith(">"):
                    seq.append(line.strip().upper())
        return "".join(seq)


class MockAlignerFactory:
    """Create multiple mock aligners with different accuracy levels."""
    
    @staticmethod
    def create_mocks(config, output_dir):
        configs = [
            ("perfect", 0.98),
            ("good", 0.88),
            ("medium", 0.75),
            ("poor", 0.55),
        ]
        return [MockAligner(config, output_dir, acc, name)
                for name, acc in configs]
