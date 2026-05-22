from .base import AlignerBase


class BasalAligner(AlignerBase):
    name = "basal"
    binary_name = "basal"
    requires_index = True

    def build_index(self, reference_fa):
        cmd = ["basal", "index", reference_fa]
        self._run_cmd(cmd, "Building basal index")

    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        out_prefix = self.aligner_dir / "aligned"
        cmd = ["basal", "--reference", reference_fa,
               "--threads", str(threads),
               read1]
        if read2:
            cmd.extend([read2])
        cmd.extend(["--output", str(out_prefix) + ".bam"])
        self._run_cmd(cmd, "Running basal alignment")
        return str(out_prefix) + ".bam"

    def call_methylation(self, sam_bam_path, reference_fa=None):
        out_path = self.aligner_dir / "methylation.bed"
        cmd = ["basal", "call", sam_bam_path, "--output", str(out_path)]
        self._run_cmd(cmd, "Calling basal methylation")
        return str(out_path)
