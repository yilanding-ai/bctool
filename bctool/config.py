import yaml
import os
from pathlib import Path
from typing import Optional


class Config:
    def __init__(self, config_path: Optional[str] = None):
        self.data = {}
        if config_path and Path(config_path).exists():
            with open(config_path) as f:
                self.data = yaml.safe_load(f)
        self._apply_defaults()

    def _apply_defaults(self):
        defaults = {
            "conversion": "ct",
            "mode": "pe",
            "input": {"read1": None, "read2": None, "reference": None},
            "output_dir": "./bctool_results",
            "extraction": {"enabled": True, "num_reads": 100000, "sample_pct": None, "seed": 42},
            "simulation": {
                "genome_length": 1000000,
                "num_reads": 50000,
                "read_length": 100,
                "fragment_mean": 300,
                "fragment_sd": 50,
                "methylation_rate": 0.7,
                "error_rate": 0.01,
                "conversion_rate": 0.99,
            },
            "aligners": {
                "enabled": [
                    "bwameth", "bsbolt", "bsmap", "walt", "abismal", "batmeth2",
                    "basal",
                    "hisat3n", "hisat3n_repeat", "bismark_bwt2_e2e", "bismark_his2",
                    "bsseeker2_bwt", "bsseeker2_soap2", "bsseeker2_bwt2_e2e",
                    "bsseeker2_bwt2_local",
                ],
                "threads": 8,
            },
            "adapter": {
                "enabled": False, "auto_detect": True,
                "r1_3prime": "AGATCGGAAGAGC",
                "r2_3prime": "AGATCGGAAGAGC",
                "r1_5prime": None, "r2_5prime": None,
                "use_cutadapt": True, "min_overlap": 6,
            },
            "umi": {
                "enabled": False, "location": "header",
                "length": 0, "pattern": None,
                "strategy": "cluster", "dedup": True,
            },
            "correction": {
                "enabled": False,
                "strategies": ["mq", "clip", "consensus"],
                "min_mq": 20,
                "max_clip_pct": 50,
                "max_unconverted": 3,
                "min_depth": 5,
                "min_agree_ratio": 0.7,
                "read_length": 150,
                "motifs": [],
                "auto_detect_motifs": True,
                "motif_flank": None,
                "prior_strength": None,
                "model_threshold": 0.5,
            },
            "evaluation": {"min_coverage": 5, "min_qual": 20},
            "threads": 8,
        }
        for key, val in defaults.items():
            if key not in self.data:
                self.data[key] = val
            elif isinstance(val, dict):
                for subkey, subval in val.items():
                    if subkey not in self.data[key]:
                        self.data[key][subkey] = subval

    def get(self, *keys, default=None):
        d = self.data
        for k in keys:
            if isinstance(d, dict) and k in d:
                d = d[k]
            else:
                return default
        return d

    @property
    def conversion(self) -> str:
        return self.get("conversion")

    @property
    def mode(self) -> str:
        return self.get("mode")

    @property
    def read1(self) -> Optional[str]:
        return self.get("input", "read1")

    @property
    def read2(self) -> Optional[str]:
        return self.get("input", "read2")

    @property
    def reference(self) -> Optional[str]:
        return self.get("input", "reference")

    @property
    def output_dir(self) -> str:
        return self.get("output_dir")

    @property
    def enabled_aligners(self) -> list:
        return self.get("aligners", "enabled", default=[])

    @property
    def threads(self) -> int:
        return self.get("threads", default=8)
