import json
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def generate_report(result_dir):
    result_dir = Path(result_dir)
    report_dir = result_dir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    results_file = result_dir / "comparison_results.json"
    if not results_file.exists():
        print(f"No results found at {results_file}")
        return

    with open(results_file) as f:
        data = json.load(f)

    df = pd.DataFrame(data["results"])
    config = data.get("config", {})

    _write_summary(report_dir, df, config)
    _plot_bar_chart(report_dir, df, "f1_score", "F1 Score by Aligner")
    _plot_bar_chart(report_dir, df, "precision", "Precision by Aligner")
    _plot_bar_chart(report_dir, df, "recall", "Recall by Aligner")
    _plot_bar_chart(report_dir, df, "accuracy", "Accuracy by Aligner")
    _plot_heatmap(report_dir, df)
    _write_html(report_dir, df, config)

    print(f"Report generated in {report_dir}")
    print(f"Open {report_dir / 'report.html'} in your browser")


def _write_summary(report_dir, df, config):
    path = report_dir / "summary.txt"
    with open(path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("BCTool Results Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Conversion type: {config.get('conversion', 'ct')}\n")
        f.write(f"Sequencing mode: {config.get('mode', 'pe')}\n")
        f.write(f"Total aligners tested: {len(df)}\n\n")

        if not df.empty:
            best = df.iloc[0]
            f.write(f"Best aligner: {best['aligner']}\n")
            f.write(f"  F1 Score:     {best['f1_score']:.4f}\n")
            f.write(f"  Precision:    {best['precision']:.4f}\n")
            f.write(f"  Recall:       {best['recall']:.4f}\n")
            f.write(f"  Accuracy:     {best['accuracy']:.4f}\n")
            f.write(f"  Level Corr:   {best.get('level_correlation', 'N/A')}\n\n")

        f.write("-" * 60 + "\n")
        f.write("Full Ranking:\n")
        f.write("-" * 60 + "\n")
        for _, row in df.iterrows():
            f.write(f"  #{row['rank']:>2}  {row['aligner']:<25} "
                   f"F1={row['f1_score']:.4f}  P={row['precision']:.4f}  "
                   f"R={row['recall']:.4f}  A={row['accuracy']:.4f}\n")

        f.write("\nDetails:\n")
        for _, row in df.iterrows():
            f.write(f"\n  {row['aligner']}:\n")
            f.write(f"    Status: {row.get('status', 'N/A')}\n")
            f.write(f"    Sites called: {row.get('sites_called', 'N/A')}\n")
            f.write(f"    Sites overlap: {row.get('sites_overlap', 'N/A')}\n")
            if 'true_positives' in row:
                f.write(f"    TP={row['true_positives']} FP={row['false_positives']} "
                       f"TN={row['true_negatives']} FN={row['false_negatives']}\n")


def _plot_bar_chart(report_dir, df, metric, title):
    if df.empty:
        return
    plt.figure(figsize=(14, 6))
    colors = sns.color_palette("coolwarm", n_colors=len(df))
    bars = plt.bar(range(len(df)), df[metric], color=colors)
    plt.xticks(range(len(df)), df["aligner"], rotation=45, ha="right", fontsize=9)
    plt.ylabel(metric.replace("_", " ").title())
    plt.title(title)
    plt.tight_layout()

    for bar, val in zip(bars, df[metric]):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                 f"{val:.4f}", ha="center", va="bottom", fontsize=8)

    plt.savefig(report_dir / f"{metric}.png", dpi=150)
    plt.close()


def _plot_heatmap(report_dir, df):
    metrics = ["f1_score", "precision", "recall", "accuracy"]
    if df.empty or not all(m in df.columns for m in metrics):
        return

    plt.figure(figsize=(10, max(4, len(df) * 0.4)))
    plot_data = df.set_index("aligner")[metrics]
    sns.heatmap(plot_data, annot=True, fmt=".4f", cmap="YlOrRd",
                linewidths=0.5, cbar_kws={"label": "Score"})
    plt.title("Aligner Performance Metrics")
    plt.tight_layout()
    plt.savefig(report_dir / "heatmap.png", dpi=150)
    plt.close()


def _write_html(report_dir, df, config):
    html = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>BCTool Report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }
  h1 { color: #333; border-bottom: 2px solid #4a90d9; padding-bottom: 10px; }
  h2 { color: #555; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; background: white;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; }
  th { background: #4a90d9; color: white; padding: 12px; text-align: left; }
  td { padding: 10px 12px; border-bottom: 1px solid #eee; }
  tr:hover { background: #f0f6ff; }
  .best { background: #e8f5e9 !important; font-weight: bold; }
  img { max-width: 100%; margin: 20px 0; box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        border-radius: 8px; }
  .summary { background: white; padding: 20px; border-radius: 8px;
             box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin: 20px 0; }
  .badge { display: inline-block; background: #4a90d9; color: white;
           padding: 4px 12px; border-radius: 12px; font-size: 14px; margin: 2px; }
  .warning { background: #fff3e0; border-left: 4px solid #ff9800; padding: 10px; }
</style>
</head>
<body>
"""
    html += f"<h1>BCTool Report</h1>"
    html += f"<div class='summary'>"
    html += f"<span class='badge'>Conversion: {config.get('conversion', 'ct')}</span> "
    html += f"<span class='badge'>Mode: {config.get('mode', 'pe')}</span> "
    html += f"<span class='badge'>Aligners: {len(df)}</span>"

    if not df.empty:
        best = df.iloc[0]
        html += f"<h2>Best Aligner: {best['aligner']}</h2>"
        html += f"<p>F1 Score: <strong>{best['f1_score']:.4f}</strong> | "
        html += f"Precision: {best['precision']:.4f} | "
        html += f"Recall: {best['recall']:.4f} | "
        html += f"Accuracy: {best['accuracy']:.4f}</p>"

    html += "</div>"

    html += "<h2>Performance Table</h2>"
    html += "<table><tr><th>Rank</th><th>Aligner</th><th>F1 Score</th>"
    html += "<th>Precision</th><th>Recall</th><th>Accuracy</th>"
    html += "<th>Level Corr</th><th>Status</th></tr>"

    for _, row in df.iterrows():
        cls = "best" if row["rank"] == 1 else ""
        html += f"<tr class='{cls}'>"
        html += f"<td>{row['rank']}</td>"
        html += f"<td>{row['aligner']}</td>"
        html += f"<td>{row['f1_score']:.4f}</td>"
        html += f"<td>{row['precision']:.4f}</td>"
        html += f"<td>{row['recall']:.4f}</td>"
        html += f"<td>{row['accuracy']:.4f}</td>"
        html += f"<td>{row.get('level_correlation', 'N/A')}</td>"
        html += f"<td>{row.get('status', 'N/A')}</td>"
        html += "</tr>"

    html += "</table>"

    for metric in ["f1_score", "precision", "recall", "accuracy"]:
        img_path = report_dir / f"{metric}.png"
        if img_path.exists():
            html += f"<h2>{metric.replace('_', ' ').title()}</h2>"
            html += f"<img src='{metric}.png' alt='{metric}'>"

    heatmap_path = report_dir / "heatmap.png"
    if heatmap_path.exists():
        html += "<h2>Heatmap</h2>"
        html += f"<img src='heatmap.png' alt='Heatmap'>"

    html += "</body></html>"

    with open(report_dir / "report.html", "w") as f:
        f.write(html)
