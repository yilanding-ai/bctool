from pathlib import Path
from .base import AlignerBase


class AbismalAligner(AlignerBase):
    name = "abismal"
    binary_name = "abismal"
    requires_index = True

    def build_index(self, reference_fa):
        out_dir = str(self.aligner_dir / "abismal_index")
        cmd = ["abismal", "-i", reference_fa, "-o", out_dir]
        self._run_cmd(cmd, "Building Abismal index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_sam = str(self.aligner_dir / "aligned.sam")
        cmd = ["abismal", "-i", reference_fa, "-r", read1,
               "-o", out_sam, "-t", str(threads)]
        if read2:
            cmd.extend(["-p", read2])
        self._run_cmd(cmd, "Running Abismal alignment")
        return out_sam

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        print(f"  [{self.name}] Methylation calling not directly supported, parsing SAM")
        return out_path
