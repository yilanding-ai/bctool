from pathlib import Path
from .base import AlignerBase


class WaltAligner(AlignerBase):
    name = "walt"
    binary_name = "walt"
    requires_index = True

    def build_index(self, reference_fa):
        cmd = ["walt", "-r", reference_fa,
               "-o", str(self.aligner_dir / "walt_index")]
        self._run_cmd(cmd, "Building Walt index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_prefix = str(self.aligner_dir / "aligned")
        cmd = ["walt", "-i", read1,
               "-o", out_prefix,
               "-r", reference_fa,
               "-t", str(threads)]
        if read2:
            cmd.extend(["-j", read2])
        self._run_cmd(cmd, "Running Walt alignment")
        return out_prefix + ".sam"

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        cmd = ["walt", "-r", reference_fa, "-i", sam_bam_path,
               "-o", str(self.aligner_dir / "methylation"),
               "-m", "call"]
        self._run_cmd(cmd, "Calling methylation with Walt")
        return out_path
