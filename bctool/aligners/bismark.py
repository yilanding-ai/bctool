from pathlib import Path
from .base import AlignerBase


class _BismarkBase(AlignerBase):
    binary_name = "bismark"

    def build_index(self, reference_fa):
        ref_dir = str(self.aligner_dir / "bismark_index")
        cmd = ["bismark_genome_preparation", ref_dir]
        self._run_cmd(cmd, "Preparing Bismark genome")

    def _common_cmd(self, read1, read2, threads, extra_args):
        out_dir = str(self.aligner_dir)
        cmd = ["bismark", "--genome", str(self.aligner_dir / "bismark_index"),
               "-o", out_dir,
               "-p", str(threads),
               "--non_directional"]
        if extra_args:
            cmd.extend(extra_args.split())
        cmd.append(read1)
        if read2:
            cmd.append(read2)
        return cmd

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        cmd = ["bismark_methylation_extractor",
               "-o", str(self.aligner_dir / "methylation"),
               "--bedGraph", sam_bam_path]
        self._run_cmd(cmd, "Extracting methylation with Bismark")
        return out_path


class BismarkBwt2E2EAligner(_BismarkBase):
    name = "bismark_bwt2_e2e"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_sam = str(self.aligner_dir / "aligned.sam")
        cmd = self._common_cmd(read1, read2, threads,
                               self.config.get("aligner_params", "bismark_bwt2_e2e", "extra_args"))
        self._run_cmd(cmd, "Running Bismark-bwt2-e2e alignment")
        bam = self.aligner_dir / f"{Path(read1).stem}_bismark_bt2_pe.bam"
        if not bam.exists():
            bam = self.aligner_dir / f"{Path(read1).stem}_bismark_bt2.bam"
        return str(bam)


class BismarkHis2Aligner(_BismarkBase):
    name = "bismark_his2"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        cmd = self._common_cmd(read1, read2, threads,
                               self.config.get("aligner_params", "bismark_his2", "extra_args"))
        cmd = [c.replace("bismark", "bismark --hisat2") for c in cmd]
        self._run_cmd(cmd, "Running Bismark-his2 alignment")
        bam = self.aligner_dir / f"{Path(read1).stem}_bismark_his2_pe.bam"
        if not bam.exists():
            bam = self.aligner_dir / f"{Path(read1).stem}_bismark_his2.bam"
        return str(bam)
