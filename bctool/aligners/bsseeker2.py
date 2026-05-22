from pathlib import Path
from .base import AlignerBase


class _BSseeker2Base(AlignerBase):
    binary_name = "bs_seeker2-align.py"

    def build_index(self, reference_fa):
        cmd = ["bs_seeker2-build.py", "-d", reference_fa,
               "-o", str(self.aligner_dir / "bsseeker2_index")]
        self._run_cmd(cmd, "Building BSseeker2 index")

    def _common_cmd(self, read1, read2, threads, aligner_flag):
        out_dir = str(self.aligner_dir)
        cmd = ["bs_seeker2-align.py", "-i", read1,
               "-d", self.config.reference,
               "-o", out_dir,
               "-p", str(threads),
               "--aligner", aligner_flag]
        if read2:
            cmd.extend(["-j", read2])
        return cmd

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = str(self.aligner_dir / "methylation.bed")
        cmd = ["bs_seeker2-call_methylation.py",
               "-i", sam_bam_path,
               "-d", reference_fa,
               "-o", str(self.aligner_dir / "methylation")]
        self._run_cmd(cmd, "Calling methylation with BSseeker2")
        return out_path


class BSseeker2BwtAligner(_BSseeker2Base):
    name = "bsseeker2_bwt"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        cmd = self._common_cmd(read1, read2, threads, "bwa")
        out_sam = str(self.aligner_dir / "aligned.sam")
        cmd.extend(["-o", out_sam])
        self._run_cmd(cmd, "Running BSseeker2-bwt alignment")
        return out_sam


class BSseeker2Soap2Aligner(_BSseeker2Base):
    name = "bsseeker2_soap2"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        cmd = self._common_cmd(read1, read2, threads, "soap2")
        out_sam = str(self.aligner_dir / "aligned.sam")
        cmd.extend(["-o", out_sam])
        self._run_cmd(cmd, "Running BSseeker2-soap2 alignment")
        return out_sam


class BSseeker2Bwt2E2EAligner(_BSseeker2Base):
    name = "bsseeker2_bwt2_e2e"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        cmd = self._common_cmd(read1, read2, threads, "bowtie2")
        out_sam = str(self.aligner_dir / "aligned.sam")
        cmd.extend(["-o", out_sam, "--bt2--end-to-end"])
        self._run_cmd(cmd, "Running BSseeker2-bwt2-e2e alignment")
        return out_sam


class BSseeker2Bwt2LocalAligner(_BSseeker2Base):
    name = "bsseeker2_bwt2_local"

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        cmd = self._common_cmd(read1, read2, threads, "bowtie2")
        out_sam = str(self.aligner_dir / "aligned.sam")
        cmd.extend(["-o", out_sam, "--bt2--local"])
        self._run_cmd(cmd, "Running BSseeker2-bwt2-local alignment")
        return out_sam
