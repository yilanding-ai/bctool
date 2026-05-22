from pathlib import Path
from .base import AlignerBase


class Batmeth2Aligner(AlignerBase):
    name = "batmeth2"
    binary_name = "batmeth2"
    requires_index = True

    def build_index(self, reference_fa):
        cmd = ["batmeth2", "index", reference_fa,
               "-o", str(self.aligner_dir / "batmeth2_index")]
        self._run_cmd(cmd, "Building Batmeth2 index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_dir = str(self.aligner_dir / "aligned")
        cmd = ["batmeth2", "-i", read1,
               "-r", reference_fa,
               "-o", out_dir,
               "-t", str(threads)]
        if read2:
            cmd.extend(["-j", read2])
        self._run_cmd(cmd, "Running Batmeth2 alignment")
        return out_dir + ".bam"

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        cmd = ["batmeth2", "call", "-i", sam_bam_path,
               "-r", reference_fa, "-o", str(self.aligner_dir / "methylation")]
        self._run_cmd(cmd, "Calling methylation with Batmeth2")
        return out_path
