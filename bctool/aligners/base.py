import subprocess
import shutil
from pathlib import Path
from abc import ABC, abstractmethod


class AlignerBase(ABC):
    name = "base"
    binary_name = ""
    requires_index = False

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.aligner_dir = self.output_dir / self.name
        self.aligner_dir.mkdir(parents=True, exist_ok=True)

    def is_available(self):
        return shutil.which(self.binary_name) is not None

    @abstractmethod
    def build_index(self, reference_fa):
        pass

    @abstractmethod
    def run_align(self, read1, read2=None, reference_fa=None, index_dir=None, threads=8):
        pass

    @abstractmethod
    def call_methylation(self, sam_bam_path, reference_fa=None):
        pass

    def _run_cmd(self, cmd, desc=""):
        print(f"  [{self.name}] {desc}")
        print(f"    {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    WARNING: command failed (code {result.returncode})")
            print(f"    stderr: {result.stderr[:500]}")
        return result
