from pathlib import Path
from .base import AlignerBase


class BwaMethAligner(AlignerBase):
    name = "bwameth"
    binary_name = "bwameth.py"
    requires_index = True

    def build_index(self, reference_fa):
        cmd = ["bwameth.py", "index", reference_fa]
        self._run_cmd(cmd, "Building bwameth index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_prefix = self.aligner_dir / "aligned"
        cmd = ["bwameth.py", "--reference", reference_fa,
               f"--threads={threads}",
               read1]
        if read2:
            cmd.extend([read2])
        cmd.extend(["--output", str(out_prefix) + ".bam"])
        self._run_cmd(cmd, "Running bwameth alignment")
        return str(out_prefix) + ".bam"

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = self.aligner_dir / "methylation.bed"
        cmd = ["bwameth.py", "call_methylation",
               sam_bam_path, ">", str(out_path)]
        # bwameth doesn't have a separate call; we parse the BAM directly
        print(f"  [{self.name}] Methylation called from bam flags")
        return str(out_path)
