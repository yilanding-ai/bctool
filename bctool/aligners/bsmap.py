from pathlib import Path
from .base import AlignerBase


class BSMAPAligner(AlignerBase):
    name = "bsmap"
    binary_name = "bsmap"
    requires_index = False

    def build_index(self, reference_fa):
        pass

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_bam = str(self.aligner_dir / "aligned.bam")
        cmd = ["bsmap", "-a", read1,
               "-d", reference_fa,
               "-o", out_bam,
               "-p", str(threads),
               "-n", "1"]
        if read2:
            cmd.extend(["-b", read2])
        self._run_cmd(cmd, "Running BSMAP alignment")
        return out_bam

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        cmd = ["bsmap", "-d", reference_fa,
               "-o", out_path,
               "-m", "5", sam_bam_path]
        self._run_cmd(cmd, "Calling methylation with BSMAP")
        return out_path
