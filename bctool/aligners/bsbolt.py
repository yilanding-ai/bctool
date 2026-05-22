from pathlib import Path
from .base import AlignerBase


class BSBoltAligner(AlignerBase):
    name = "bsbolt"
    binary_name = "bsbolt"
    requires_index = True

    def build_index(self, reference_fa):
        cmd = ["bsbolt", "Index", "-G", reference_fa,
               "-D", str(self.aligner_dir / "bsbolt_index")]
        self._run_cmd(cmd, "Building BSBolt index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_prefix = str(self.aligner_dir / "aligned")
        cmd = ["bsbolt", "Align", "-db", str(self.aligner_dir / "bsbolt_index"),
               "-t", str(threads), "-F1", read1,
               "-O", out_prefix]
        if read2:
            cmd.extend(["-F2", read2])
        self._run_cmd(cmd, "Running BSBolt alignment")
        return out_prefix + ".bam"

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = self.aligner_dir / "methylation.bed"
        cmd = ["bsbolt", "CallMethylation", "-I", sam_bam_path,
               "-O", str(self.aligner_dir / "methylation"),
               "-t", str(self.config.threads)]
        self._run_cmd(cmd, "Calling methylation with BSBolt")
        return str(out_path)
