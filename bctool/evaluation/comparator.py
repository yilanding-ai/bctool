import pandas as pd
from pathlib import Path


def compare_aligners(results_list):
    """Compare results from all aligners and rank them."""
    df = pd.DataFrame(results_list)

    if df.empty:
        return df

    # Sort by F1 score descending
    df = df.sort_values("f1_score", ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)

    return df


def find_best_aligner(df):
    """Return the best aligner info."""
    if df.empty:
        return None
    best = df.iloc[0]
    return {
        "aligner": best["aligner"],
        "f1_score": best["f1_score"],
        "precision": best["precision"],
        "recall": best["recall"],
        "accuracy": best["accuracy"],
        "level_correlation": best.get("level_correlation", "N/A"),
    }
