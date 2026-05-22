import pandas as pd
import numpy as np
from pathlib import Path


class Evaluator:
    """Compare aligner methylation calls against ground truth."""

    def __init__(self, ground_truth_bed):
        self.ground_truth = self._load_bed(ground_truth_bed)

    @staticmethod
    def _load_bed(path):
        data = pd.read_csv(path, sep="\t", header=None,
                           names=["chrom", "start", "end", "methylated"])
        data["pos"] = data["start"]
        data["methylated"] = data["methylated"].astype(int)
        return data

    def evaluate(self, called_bed_path, aligner_name):
        if not Path(called_bed_path).exists():
            return {
                "aligner": aligner_name,
                "precision": 0.0,
                "recall": 0.0,
                "f1_score": 0.0,
                "accuracy": 0.0,
                "sites_called": 0,
                "sites_truth": len(self.ground_truth),
                "status": "no_output"
            }

        called = pd.read_csv(called_bed_path, sep="\t", header=None,
                             names=["chrom", "start", "end", "level",
                                    "meth_count", "unmeth_count", "coverage"])
        called["pos"] = called["start"]
        called["called_meth"] = (called["level"] > 0.5).astype(int)

        merged = self.ground_truth.merge(
            called, on=["chrom", "pos"], how="inner", suffixes=("_truth", "_called"))

        if len(merged) == 0:
            return {
                "aligner": aligner_name,
                "precision": 0.0,
                "recall": 0.0,
                "f1_score": 0.0,
                "accuracy": 0.0,
                "sites_called": 0,
                "sites_truth": len(self.ground_truth),
                "sites_overlap": 0,
                "status": "no_overlap"
            }

        # TP: called methylated & truth methylated
        tp = ((merged["called_meth"] == 1) & (merged["methylated"] == 1)).sum()
        fp = ((merged["called_meth"] == 1) & (merged["methylated"] == 0)).sum()
        tn = ((merged["called_meth"] == 0) & (merged["methylated"] == 0)).sum()
        fn = ((merged["called_meth"] == 0) & (merged["methylated"] == 1)).sum()

        # methylation level correlation
        level_corr = merged["level"].corr(merged["methylated"])

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        accuracy = (tp + tn) / (tp + fp + tn + fn) if (tp + fp + tn + fn) > 0 else 0.0

        return {
            "aligner": aligner_name,
            "precision": round(precision, 4),
            "recall": round(recall, 4),
            "f1_score": round(f1, 4),
            "accuracy": round(accuracy, 4),
            "level_correlation": round(level_corr, 4),
            "true_positives": int(tp),
            "false_positives": int(fp),
            "true_negatives": int(tn),
            "false_negatives": int(fn),
            "sites_called": len(called),
            "sites_truth": len(self.ground_truth),
            "sites_overlap": len(merged),
            "status": "ok"
        }
