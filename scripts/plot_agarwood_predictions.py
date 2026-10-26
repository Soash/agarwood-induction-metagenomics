#!/usr/bin/env python3
"""
plot_agarwood_predictions.py

Visualizes agarwood_predictions_2.csv: for each fungal species, its total
Glycoside Hydrolase (GH) count, total Auxiliary Activity (AA) count, and a
predicted inducer class (Brute Force Inducer / Wound Opportunist Inducer /
Non-Inducer (Rejection)).

Since the dataset is fundamentally "two numeric features + one predicted
category, per species", the charts here are chosen to show:
  (a) how many species fall into each predicted class,
  (b) how GH and AA together separate the three classes (the two features
      driving/associated with the prediction),
  (c) how each feature's distribution differs across classes,
  (d) which specific species are the most extreme enzyme producers.

Outputs (written to OUTPUT_DIR):
    class_distribution_barchart.png     species count per predicted class
    gh_vs_aa_scatter.png                GH vs AA, colored by class (main plot)
    gh_vs_aa_scatter_marginal.png       same scatter, with marginal histograms
    gh_aa_boxplots_by_class.png         GH and AA distributions per class (2-panel)
    top15_species_by_total_enzymes.png  top 15 species by GH+AA, colored by class
    class_summary_stats.csv             mean/median/std GH & AA per class

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/plot_agarwood_predictions.py
"""

import csv
import statistics
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
except ImportError:
    plt = None

# ============================================================
# CONFIG
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/
INPUT_CSV = PROJECT_ROOT / "awi_pred" / "agarwood_predictions_2.csv"
OUTPUT_DIR = PROJECT_ROOT / "awi_pred" / "agarwood_predictions_plots"

CLASS_DIST_OUT = OUTPUT_DIR / "class_distribution_barchart.png"
SCATTER_OUT = OUTPUT_DIR / "gh_vs_aa_scatter.png"
SCATTER_MARGINAL_OUT = OUTPUT_DIR / "gh_vs_aa_scatter_marginal.png"
BOXPLOTS_OUT = OUTPUT_DIR / "gh_aa_boxplots_by_class.png"
TOP15_OUT = OUTPUT_DIR / "top15_species_by_total_enzymes.png"
SUMMARY_CSV_OUT = OUTPUT_DIR / "class_summary_stats.csv"

CLASS_ORDER = ["Brute Force Inducer", "Wound Opportunist Inducer", "Non-Inducer (Rejection)"]
CLASS_COLORS = {
    "Brute Force Inducer": "#C44E52",
    "Wound Opportunist Inducer": "#DD8452",
    "Non-Inducer (Rejection)": "#8C8C8C",
}
TOP_N = 15


# ============================================================
# STEP 1: Load the data
# ============================================================

def load_data(path: Path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                gh = int(row["Total GH"])
                aa = int(row["Total AA"])
            except (ValueError, KeyError):
                continue
            rows.append({
                "species": row["Species"],
                "gh": gh,
                "aa": aa,
                "prediction": row["Prediction"].strip(),
            })
    return rows


# ============================================================
# STEP 2: Class distribution bar chart
# ============================================================

def plot_class_distribution(rows):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping class distribution chart.")
        return

    counts = {c: 0 for c in CLASS_ORDER}
    for r in rows:
        counts[r["prediction"]] = counts.get(r["prediction"], 0) + 1

    labels = [c for c in CLASS_ORDER if c in counts]
    values = [counts[c] for c in labels]
    colors = [CLASS_COLORS.get(c, "#4C72B0") for c in labels]

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(labels, values, color=colors)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(val),
                ha="center", va="bottom")

    total = sum(values)
    ax.set_ylabel("Number of species")
    ax.set_title(f"Predicted Inducer Class Distribution (n={total} species)")
    plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
    fig.savefig(CLASS_DIST_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Class distribution bar chart saved to: {CLASS_DIST_OUT}")


# ============================================================
# STEP 3a: GH vs AA scatter, colored by class
# ============================================================

def plot_scatter(rows):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping scatter plot.")
        return

    fig, ax = plt.subplots(figsize=(9, 7))
    for cls in CLASS_ORDER:
        xs = [r["gh"] for r in rows if r["prediction"] == cls]
        ys = [r["aa"] for r in rows if r["prediction"] == cls]
        ax.scatter(xs, ys, label=f"{cls} (n={len(xs)})",
                   color=CLASS_COLORS.get(cls), alpha=0.65, edgecolors="none", s=35)

    ax.set_xlabel("Total GH (Glycoside Hydrolase) count")
    ax.set_ylabel("Total AA (Auxiliary Activity) count")
    ax.set_title("Fungal Species: GH vs AA Enzyme Counts by Predicted Class")
    ax.legend(loc="upper left")
    fig.savefig(SCATTER_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  GH vs AA scatter plot saved to: {SCATTER_OUT}")


# ============================================================
# STEP 3b: Same scatter, with marginal histograms
# ============================================================

def plot_scatter_with_marginals(rows):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping marginal scatter plot.")
        return

    fig = plt.figure(figsize=(9, 9))
    gs = GridSpec(4, 4, figure=fig, hspace=0.05, wspace=0.05)
    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    for cls in CLASS_ORDER:
        xs = [r["gh"] for r in rows if r["prediction"] == cls]
        ys = [r["aa"] for r in rows if r["prediction"] == cls]
        color = CLASS_COLORS.get(cls)
        ax_main.scatter(xs, ys, label=cls, color=color, alpha=0.6, edgecolors="none", s=30)
        ax_top.hist(xs, bins=25, color=color, alpha=0.5)
        ax_right.hist(ys, bins=25, color=color, alpha=0.5, orientation="horizontal")

    ax_main.set_xlabel("Total GH count")
    ax_main.set_ylabel("Total AA count")
    ax_main.legend(loc="upper left", fontsize=8)

    plt.setp(ax_top.get_xticklabels(), visible=False)
    plt.setp(ax_right.get_yticklabels(), visible=False)
    ax_top.set_ylabel("Count")
    ax_right.set_xlabel("Count")
    ax_top.set_title("GH vs AA by Predicted Class (with marginal distributions)")

    fig.savefig(SCATTER_MARGINAL_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Scatter with marginal histograms saved to: {SCATTER_MARGINAL_OUT}")


# ============================================================
# STEP 4: Box plots of GH and AA by class (2-panel)
# ============================================================

def plot_boxplots(rows):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping box plots.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for ax, feature, ylabel in [(axes[0], "gh", "Total GH count"), (axes[1], "aa", "Total AA count")]:
        data = [[r[feature] for r in rows if r["prediction"] == cls] for cls in CLASS_ORDER]
        bp = ax.boxplot(data, tick_labels=[c.replace(" (", "\n(") for c in CLASS_ORDER], patch_artist=True)
        for patch, cls in zip(bp["boxes"], CLASS_ORDER):
            patch.set_facecolor(CLASS_COLORS.get(cls))
            patch.set_alpha(0.7)
        ax.set_ylabel(ylabel)
        plt.setp(ax.get_xticklabels(), rotation=0, fontsize=8)

    fig.suptitle("Enzyme Count Distributions by Predicted Class")
    fig.tight_layout()
    fig.savefig(BOXPLOTS_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Box plots saved to: {BOXPLOTS_OUT}")


# ============================================================
# STEP 5: Top N species by total enzymes (GH+AA), colored by class
# ============================================================

def plot_top_species(rows):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping top-species chart.")
        return

    ranked = sorted(rows, key=lambda r: -(r["gh"] + r["aa"]))[:TOP_N]
    ranked = ranked[::-1]  # so the highest ends up at the top of a horizontal bar chart

    species = [r["species"] for r in ranked]
    totals = [r["gh"] + r["aa"] for r in ranked]
    colors = [CLASS_COLORS.get(r["prediction"], "#4C72B0") for r in ranked]

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.barh(species, totals, color=colors)
    ax.set_xlabel("Total enzyme count (GH + AA)")
    ax.set_title(f"Top {TOP_N} Species by Total Enzyme Count")

    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c]) for c in CLASS_ORDER]
    ax.legend(handles, CLASS_ORDER, loc="lower right", fontsize=8)

    fig.savefig(TOP15_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Top {TOP_N} species bar chart saved to: {TOP15_OUT}")


# ============================================================
# STEP 6: Summary stats CSV per class
# ============================================================

def write_summary_csv(rows):
    with open(SUMMARY_CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Class", "n_species", "mean_GH", "median_GH", "stdev_GH",
                          "mean_AA", "median_AA", "stdev_AA"])
        for cls in CLASS_ORDER:
            gh_vals = [r["gh"] for r in rows if r["prediction"] == cls]
            aa_vals = [r["aa"] for r in rows if r["prediction"] == cls]
            if not gh_vals:
                continue
            writer.writerow([
                cls, len(gh_vals),
                round(statistics.mean(gh_vals), 2), statistics.median(gh_vals),
                round(statistics.stdev(gh_vals), 2) if len(gh_vals) > 1 else 0,
                round(statistics.mean(aa_vals), 2), statistics.median(aa_vals),
                round(statistics.stdev(aa_vals), 2) if len(aa_vals) > 1 else 0,
            ])
    print(f"  Class summary stats CSV saved to: {SUMMARY_CSV_OUT}")


# ============================================================
# MAIN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Could not find {INPUT_CSV}. Update INPUT_CSV at the top of this script.")

    rows = load_data(INPUT_CSV)
    print(f"Loaded {len(rows)} species from {INPUT_CSV.name}\n")

    print("=" * 60)
    print("Generating charts")
    print("=" * 60)
    plot_class_distribution(rows)
    plot_scatter(rows)
    plot_scatter_with_marginals(rows)
    plot_boxplots(rows)
    plot_top_species(rows)
    write_summary_csv(rows)


if __name__ == "__main__":
    main()
    
    
    