"""
Multi-sample batch analysis module.

Enables running the BCTool pipeline across multiple samples
with different conditions, then comparing results to find which
aligner performs best under each condition.
"""

import json
import yaml
import copy
import time
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

import pandas as pd

from .config import Config
from .methylation.converter import ConversionType
from .pipeline import Pipeline
from .aligners.mock import MockAlignerFactory


class MultiSampleRunner:
    """Run BCTool across multiple samples and compare results."""

    def __init__(self, config_path=None):
        self.config_path = config_path
        self.samples = []
        self.results = {}
        self._load_config()

    def _load_config(self):
        if self.config_path and Path(self.config_path).exists():
            with open(self.config_path) as f:
                data = yaml.safe_load(f)
            self.samples = data.get("samples", [])
            self.global_config = {k: v for k, v in data.items() if k != "samples"}
        else:
            self.samples = []
            self.global_config = {}

    def add_sample(self, name, read1, read2=None, reference=None,
                   conversion="ct", mode="pe", group=None,
                   sim_reads=0, **kwargs):
        """Add a sample to the batch."""
        sample = {
            "name": name,
            "read1": read1,
            "read2": read2,
            "reference": reference,
            "conversion": conversion,
            "mode": mode,
            "group": group or "default",
            "sim_reads": sim_reads,
        }
        sample.update(kwargs)
        self.samples.append(sample)

    def run_all(self, output_dir, max_workers=2, use_mock=False):
        """Run pipeline on all samples, optionally in parallel."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        start = datetime.now()

        print("=" * 70)
        print(f"  BCTool Multi-Sample Analysis")
        print(f"  Samples: {len(self.samples)}  Workers: {max_workers}")
        print(f"  Started: {start}")
        print("=" * 70)

        def _run_sample(sample):
            name = sample["name"]
            print(f"\n>>> Processing sample: {name} <<<")
            sample_out = out / name
            sample_out.mkdir(parents=True, exist_ok=True)

            cfg = Config()
            for k, v in self.global_config.items():
                if isinstance(v, dict):
                    for sk, sv in v.items():
                        cfg.data.setdefault(k, {})[sk] = sv
                else:
                    cfg.data[k] = v

            cfg.data["output_dir"] = str(sample_out)
            cfg.data["conversion"] = sample.get("conversion", "ct")
            cfg.data["mode"] = sample.get("mode", "pe")
            cfg.data["simulation"]["num_reads"] = sample.get("sim_reads", 0)
            cfg.data["extraction"]["enabled"] = False

            if sample.get("sim_reads", 0) > 0:
                pass  # will simulate
            else:
                cfg.data["input"]["read1"] = sample["read1"]
                if sample.get("read2"):
                    cfg.data["input"]["read2"] = sample["read2"]
                if sample.get("reference"):
                    cfg.data["input"]["reference"] = sample["reference"]

            if use_mock:
                cfg.data["aligners"]["enabled"] = []

            try:
                pipe = Pipeline(cfg)
                pipe.run()
                res_file = sample_out / "comparison_results.json"
                if res_file.exists():
                    with open(res_file) as f:
                        result = json.load(f)
                    return name, result, None
                return name, None, "no results file"
            except Exception as e:
                return name, None, str(e)

        self.results = {}
        errors = {}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_run_sample, s): s["name"]
                       for s in self.samples}
            for future in as_completed(futures):
                name, result, error = future.result()
                if error:
                    errors[name] = error
                    print(f"  [ERROR] {name}: {error}")
                else:
                    self.results[name] = result

        # Save raw results
        with open(out / "multi_sample_results.json", "w") as f:
            json.dump({
                "samples": self.samples,
                "results": {k: v for k, v in self.results.items() if v},
                "errors": errors,
                "timestamp": str(start),
            }, f, indent=2)

        elapsed = (datetime.now() - start).total_seconds()
        print(f"\n{'=' * 70}")
        print(f"  Multi-sample analysis complete")
        print(f"  Successful: {len(self.results)}/{len(self.samples)}")
        print(f"  Errors: {len(errors)}")
        print(f"  Time: {elapsed:.1f}s")
        print(f"{'=' * 70}")

        # Generate comparison
        if self.results:
            self._compare_results(out)

        return self.results, errors

    def _compare_results(self, output_dir):
        """Cross-sample comparison and ranking."""
        print("\n>>> Generating cross-sample comparison...")

        # Build comparison table
        rows = []
        for sample_name, result in self.results.items():
            if not result or "results" not in result:
                continue
            sample_info = next((s for s in self.samples if s["name"] == sample_name), {})
            group = sample_info.get("group", "default")
            conversion = sample_info.get("conversion", result.get("config", {}).get("conversion", "ct"))

            for r in result["results"]:
                rows.append({
                    "sample": sample_name,
                    "group": group,
                    "conversion": conversion,
                    "aligner": r["aligner"],
                    "f1_score": r.get("f1_score", 0),
                    "precision": r.get("precision", 0),
                    "recall": r.get("recall", 0),
                    "accuracy": r.get("accuracy", 0),
                    "sites_overlap": r.get("sites_overlap", 0),
                    "status": r.get("status", ""),
                })

        if not rows:
            print("  No results to compare")
            return

        df = pd.DataFrame(rows)

        # 1. Overall winner per sample
        print("\n  --- Best Aligner Per Sample ---")
        best_per_sample = df.loc[df.groupby("sample")["f1_score"].idxmax()]
        for _, row in best_per_sample.iterrows():
            print(f"    {row['sample']:<20} -> {row['aligner']:<25} "
                  f"(F1={row['f1_score']:.4f})")

        # 2. Overall winner across all samples
        print("\n  --- Overall Aligner Ranking (all samples) ---")
        overall = df.groupby("aligner").agg(
            mean_f1=("f1_score", "mean"),
            std_f1=("f1_score", "std"),
            wins=("f1_score", lambda x: (x == x.max()).sum() if len(x) > 0 else 0),
            count=("f1_score", "count"),
        ).sort_values("mean_f1", ascending=False)
        overall["rank"] = range(1, len(overall) + 1)
        print(f"  {'Rank':<6} {'Aligner':<25} {'Mean F1':<10} {'Std':<8} {'Wins':<6}")
        print("  " + "-" * 55)
        for _, row in overall.iterrows():
            print(f"  {int(row['rank']):<6} {row.name:<25} "
                  f"{row['mean_f1']:<10.4f} {row['std_f1']:<8.4f} {int(row['wins']):<6}")

        # 3. Per-group analysis
        if df["group"].nunique() > 1:
            print("\n  --- Per-Group Best Aligner ---")
            for group in df["group"].unique():
                gdf = df[df["group"] == group]
                grank = gdf.groupby("aligner")["f1_score"].mean().sort_values(ascending=False)
                print(f"\n    Group: {group}")
                for i, (al, f1) in enumerate(grank.items(), 1):
                    marker = " <<< BEST" if i == 1 else ""
                    print(f"      #{i:<3} {al:<25} F1={f1:.4f}{marker}")

        # 4. Per-conversion analysis
        if df["conversion"].nunique() > 1:
            print("\n  --- Per-Conversion Best Aligner ---")
            for conv in df["conversion"].unique():
                cdf = df[df["conversion"] == conv]
                crank = cdf.groupby("aligner")["f1_score"].mean().sort_values(ascending=False)
                conv_label = ConversionType(conv).key
                print(f"\n    Conversion: {conv_label} ({ConversionType(conv).description()})")
                for i, (al, f1) in enumerate(crank.items(), 1):
                    marker = " <<< BEST" if i == 1 else ""
                    print(f"      #{i:<3} {al:<25} F1={f1:.4f}{marker}")

        # Save comparison tables
        df.to_csv(output_dir / "multi_sample_comparison.csv", index=False)
        overall.to_csv(output_dir / "overall_ranking.csv")

        # Summary
        summary = {
            "total_samples": len(self.results),
            "total_aligners": df["aligner"].nunique(),
            "overall_best": overall.index[0] if len(overall) > 0 else None,
            "overall_best_f1": float(overall.iloc[0]["mean_f1"]) if len(overall) > 0 else 0,
            "per_sample_best": best_per_sample[["sample", "aligner", "f1_score"]].to_dict("records"),
        }
        with open(output_dir / "multi_sample_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        # Generate comparison report
        self._generate_comparison_report(output_dir, df, overall, best_per_sample)
        return df, overall

    def _generate_comparison_report(self, output_dir, df, overall, best_per_sample):
        """Generate HTML report for multi-sample comparison."""
        report_dir = Path(output_dir) / "comparison_report"
        report_dir.mkdir(exist_ok=True)

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import seaborn as sns
            has_plot = True
        except ImportError:
            has_plot = False

        if has_plot:
            # 1. Per-sample bar chart
            plt.figure(figsize=(14, 8))
            pivot = df.pivot_table(index="sample", columns="aligner",
                                    values="f1_score", aggfunc="mean")
            pivot.plot(kind="bar", ax=plt.gca())
            plt.title("F1 Score by Aligner per Sample", fontsize=14)
            plt.ylabel("F1 Score")
            plt.xlabel("Sample")
            plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(report_dir / "per_sample_f1.png", dpi=150)
            plt.close()

            # 2. Heatmap
            plt.figure(figsize=(10, max(5, len(df["sample"].unique()) * 0.4)))
            heat = df.pivot_table(index="sample", columns="aligner",
                                   values="f1_score", aggfunc="mean")
            sns.heatmap(heat, annot=True, fmt=".3f", cmap="YlOrRd",
                        linewidths=0.5, cbar_kws={"label": "F1 Score"})
            plt.title("F1 Score Heatmap (Sample x Aligner)", fontsize=14)
            plt.tight_layout()
            plt.savefig(report_dir / "comparison_heatmap.png", dpi=150)
            plt.close()

            # 3. Win count chart
            if "wins" in overall.columns:
                plt.figure(figsize=(10, 6))
                colors = sns.color_palette("coolwarm", n_colors=len(overall))
                plt.bar(range(len(overall)), overall["wins"], color=colors)
                plt.xticks(range(len(overall)), overall.index, rotation=45, ha="right")
                plt.title("Number of Times Each Aligner Won (Best F1)")
                plt.ylabel("Win Count")
                for i, w in enumerate(overall["wins"]):
                    plt.text(i, w + 0.1, str(int(w)), ha="center")
                plt.tight_layout()
                plt.savefig(report_dir / "win_counts.png", dpi=150)
                plt.close()

        # HTML report
        html = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Multi-Sample Comparison Report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }
  h1 { color: #333; border-bottom: 2px solid #4a90d9; padding-bottom: 10px; }
  h2 { color: #555; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; background: white;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; }
  th { background: #4a90d9; color: white; padding: 10px; text-align: left; }
  td { padding: 8px 10px; border-bottom: 1px solid #eee; }
  tr:hover { background: #f0f6ff; }
  .best { background: #e8f5e9 !important; font-weight: bold; }
  img { max-width: 100%; margin: 15px 0; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
  .summary-box { background: white; padding: 20px; border-radius: 8px;
                 box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin: 20px 0; }
  .badge { display: inline-block; padding: 4px 12px; border-radius: 12px;
           font-size: 13px; margin: 2px; color: white; }
  .badge-blue { background: #4a90d9; }
  .badge-green { background: #43a047; }
</style>
</head>
<body>
"""
        html += "<h1>Multi-Sample Comparison Report</h1>"

        html += "<div class='summary-box'>"
        html += f"<span class='badge badge-blue'>Samples: {df['sample'].nunique()}</span> "
        html += f"<span class='badge badge-blue'>Aligners: {df['aligner'].nunique()}</span> "
        if "overall_best" in df.columns:
            pass
        html += "</div>"

        # Overall ranking table
        html += "<h2>Overall Aligner Ranking</h2>"
        html += "<table><tr><th>Rank</th><th>Aligner</th><th>Mean F1</th><th>Std F1</th><th>Wins</th></tr>"
        for _, row in overall.iterrows():
            cls = "best" if row["rank"] == 1 else ""
            html += f"<tr class='{cls}'>"
            html += f"<td>{int(row['rank'])}</td>"
            html += f"<td>{row.name}</td>"
            html += f"<td>{row['mean_f1']:.4f}</td>"
            html += f"<td>{row['std_f1']:.4f}</td>"
            html += f"<td>{int(row['wins'])}</td>"
            html += "</tr>"
        html += "</table>"

        # Per-sample best
        html += "<h2>Best Aligner per Sample</h2>"
        html += "<table><tr><th>Sample</th><th>Best Aligner</th><th>F1 Score</th></tr>"
        for _, row in best_per_sample.iterrows():
            html += f"<tr><td>{row['sample']}</td><td><strong>{row['aligner']}</strong></td>"
            html += f"<td>{row['f1_score']:.4f}</td></tr>"
        html += "</table>"

        # Per-conversion analysis
        if "conversion" in df.columns and df["conversion"].nunique() > 1:
            html += "<h2>Per-Conversion Analysis</h2>"
            for conv in df["conversion"].unique():
                cdf = df[df["conversion"] == conv]
                crank = cdf.groupby("aligner")["f1_score"].mean().sort_values(ascending=False)
                label = ConversionType(conv).key
                html += f"<h3>Conversion: {label} ({ConversionType(conv).description()})</h3>"
                html += "<table><tr><th>Rank</th><th>Aligner</th><th>Mean F1</th></tr>"
                for i, (al, f1) in enumerate(crank.items(), 1):
                    cls = "best" if i == 1 else ""
                    html += f"<tr class='{cls}'><td>{i}</td><td>{al}</td><td>{f1:.4f}</td></tr>"
                html += "</table>"

        # Charts
        for img in ["per_sample_f1.png", "comparison_heatmap.png", "win_counts.png"]:
            img_path = report_dir / img
            if img_path.exists():
                html += f"<h2>{img.replace('_', ' ').replace('.png', '').title()}</h2>"
                html += f"<img src='{img}' alt='{img}'>"

        html += "</body></html>"

        with open(report_dir / "comparison_report.html", "w") as f:
            f.write(html)

        print(f"  Comparison report: {report_dir / 'comparison_report.html'}")


def generate_batch_config(output_path, samples=None):
    """Generate a sample multi-sample batch config file."""
    if samples is None:
        samples = [
            {
                "name": "sample_ct_high_meth",
                "conversion": "ct",
                "mode": "pe",
                "group": "high_methylation",
                "sim_reads": 30000,
            },
            {
                "name": "sample_ct_low_meth",
                "conversion": "ct",
                "mode": "pe",
                "group": "low_methylation",
                "sim_reads": 30000,
            },
            {
                "name": "sample_ag_high_meth",
                "conversion": "ag",
                "mode": "pe",
                "group": "high_methylation",
                "sim_reads": 30000,
            },
            {
                "name": "sample_ag_low_meth",
                "conversion": "ag",
                "mode": "pe",
                "group": "low_methylation",
                "sim_reads": 30000,
            },
        ]

    config = {
        "max_workers": 2,
        "use_mock": True,
        "samples": samples,
    }

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"Batch config written to: {path}")
    return str(path)
