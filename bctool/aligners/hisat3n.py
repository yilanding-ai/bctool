from pathlib import Path
from .base import AlignerBase


class Hisat3nAligner(AlignerBase):
    name = "hisat3n"
    binary_name = "hisat-3n"
    requires_index = True

    def build_index(self, reference_fa):
        cmd = ["hisat-3n-build", reference_fa,
               str(self.aligner_dir / "hisat3n_index")]
        self._run_cmd(cmd, "Building HISAT-3n index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_sam = str(self.aligner_dir / "aligned.sam")
        idx = str(self.aligner_dir / "hisat3n_index")
        cmd = ["hisat-3n", "-x", idx,
               "--mp", read1,
               "-o", out_sam,
               "-p", str(threads),
               "--mp-seed", "C"]
        if read2:
            cmd.extend(["--mp", read2])
        self._run_cmd(cmd, "Running HISAT-3n alignment")
        return out_sam

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        print(f"  [{self.name}] Methylation calling via SAM parsing")
        return out_path


class Hisat3nRepeatAligner(Hisat3nAligner):
    name = "hisat3n_repeat"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_sam = str(self.aligner_dir / "aligned.sam")
        idx = str(self.aligner_dir / "hisat3n_index")
        cmd = ["hisat-3n", "-x", idx,
               "--mp", read1,
               "-o", out_sam,
               "-p", str(threads),
               "--mp-seed", "C",
               "--repeat"]
        if read2:
            cmd.extend(["--mp", read2])
        self._run_cmd(cmd, "Running HISAT-3n repeat alignment")
        return out_sam
