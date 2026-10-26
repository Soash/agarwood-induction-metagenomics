#!/usr/bin/env python3
"""
fungi_species_heatmap.py

Clustered heatmap of fungal species abundance (CPM) across samples.

Rows = top N fungal species (ranked by the same weighted_abundance_score
used in fungi_summary.py: balances how widespread a species is across the
non-control samples with how abundant it is when present).
Columns = all samples, including Control (kept in for visual contrast --
a near-empty Control column is itself informative).
Cell color = log10(CPM + 1), since CPM values span several orders of
magnitude; a linear color scale would make everything except the single
brightest cell look the same.

Rows and columns are hierarchically clustered (average-linkage, Euclidean
distance on the log-CPM matrix) so species with similar abundance patterns
across samples sit next to each other, and samples with similar fungal
communities sit next to each other -- this is what turns a plain grid into
a genuinely informative "clustermap".

Uses the same Fungi-subtree parsing and EXCLUDED_SPECIES contaminant list
as fungi_summary.py, so results stay consistent across scripts.

Outputs (written to OUTPUT_DIR):
    fungi_species_heatmap.png   the clustered heatmap with dendrograms

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/fungi_species_heatmap.py
"""

import math
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from scipy.cluster.hierarchy import linkage, dendrogram
except ImportError:
    plt = None

# ============================================================
# CONFIG — adjust these if your layout differs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/

BRACKEN_REPORTS_DIR = PROJECT_ROOT / "bracken_results"
OUTPUT_DIR = PROJECT_ROOT / "fungi_summary"

HEATMAP_OUT = OUTPUT_DIR / "fungi_species_heatmap.png"

TAXID_FUNGI = 4751
TAXID_ROOT = 1

CONTROL_SAMPLE = "Control"
TOP_N = 20  # number of species shown as rows

PREVALENCE_WEIGHT = 0.5  # same as fungi_summary.py's weighted_abundance_score

EXCLUDED_SPECIES = {
    "Malassezia globosa",
    "Malassezia restricta",
    "Malassezia japonica",
    "uncultured fungus",
    "uncultured Ascomycota",
    "fungal sp.",
}

# Samples are auto-detected from bracken_results/*_bracken_report.txt,
# but you can hardcode instead if you prefer:
# SAMPLES = ["AI_1", "AI_2", "AI_3", "KNI", "NI_1", "NI_2", "Control"]
SAMPLES = None  # None = auto-detect


# ============================================================
# STEP 1: Parse one report -> (root_total, {species_name: count})
# (identical logic to fungi_summary.py's parse_fungi_species)
# ============================================================

def parse_fungi_species(path: Path):
    root_total = 0
    species_counts = {}

    within_fungi = False
    fungi_indent = None

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue

            try:
                count_clade = int(fields[1])
                rank_code = fields[3]
                tax_id = int(fields[4])
            except ValueError:
                continue

            name_field = fields[5]
            indent = len(name_field) - len(name_field.lstrip(" "))
            name = name_field.strip()

            if tax_id == TAXID_ROOT:
                root_total = count_clade

            if not within_fungi:
                if tax_id == TAXID_FUNGI:
                    within_fungi = True
                    fungi_indent = indent
                continue

            if indent <= fungi_indent:
                within_fungi = False
                continue

            if rank_code == "S":
                if name in EXCLUDED_SPECIES:
                    continue
                species_counts[name] = species_counts.get(name, 0) + count_clade

    return root_total, species_counts


# ============================================================
# STEP 2: Rank species by weighted_abundance_score
# (identical logic to fungi_summary.py's rank_species_by_weighted_score)
# ============================================================

def rank_species_by_weighted_score(all_species, all_root_totals, samples, top_n):
    non_control_samples = [s for s in samples if s != CONTROL_SAMPLE]

    all_species_names = set()
    for sample in non_control_samples:
        all_species_names.update(all_species.get(sample, {}).keys())

    prevalence = {}
    mean_cpm = {}
    for species in all_species_names:
        cpm_vals = []
        for sample in non_control_samples:
            raw = all_species.get(sample, {}).get(species, 0)
            root = all_root_totals.get(sample, 0)
            if raw > 0 and root:
                cpm_vals.append(raw / root * 1_000_000)
        prevalence[species] = len(cpm_vals)
        mean_cpm[species] = sum(cpm_vals) / len(cpm_vals) if cpm_vals else 0.0

    max_mean = max(mean_cpm.values(), default=0.0)
    scores = {}
    for species in all_species_names:
        prevalence_fraction = prevalence[species] / len(non_control_samples) if non_control_samples else 0.0
        abundance_fraction = mean_cpm[species] / max_mean if max_mean > 0 else 0.0
        scores[species] = (PREVALENCE_WEIGHT * prevalence_fraction
                            + (1 - PREVALENCE_WEIGHT) * abundance_fraction)

    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    return [s for s, _ in ranked[:top_n]]


# ============================================================
# STEP 3: Build the CPM matrix (log10(CPM+1)) for the top species
# ============================================================

def build_matrix(top_species, samples, all_species, all_root_totals):
    matrix = []
    for species in top_species:
        row = []
        for sample in samples:
            raw = all_species.get(sample, {}).get(species, 0)
            root = all_root_totals.get(sample, 0)
            cpm = (raw / root * 1_000_000) if root else 0.0
            row.append(math.log10(cpm + 1))
        matrix.append(row)
    return matrix


# ============================================================
# STEP 4: Draw the clustered heatmap
# ============================================================

def plot_clustered_heatmap(matrix, top_species, samples):
    if plt is None:
        print("  WARNING: matplotlib/scipy not installed, skipping heatmap.")
        print("  Install with: pip install matplotlib scipy")
        return

    n_rows = len(matrix)
    n_cols = len(matrix[0]) if matrix else 0

    # --- Hierarchical clustering ---
    row_linkage = linkage(matrix, method="average", metric="euclidean")
    col_linkage = linkage(list(zip(*matrix)), method="average", metric="euclidean")

    row_dendro = dendrogram(row_linkage, no_plot=True)
    col_dendro = dendrogram(col_linkage, no_plot=True)
    row_order = row_dendro["leaves"]
    col_order = col_dendro["leaves"]

    ordered_matrix = [[matrix[r][c] for c in col_order] for r in row_order]
    ordered_species = [top_species[r] for r in row_order]
    ordered_samples = [samples[c] for c in col_order]

    # --- Figure layout: row dendrogram | heatmap+colorbar, with col dendrogram on top ---
    fig = plt.figure(figsize=(10, 0.35 * n_rows + 3))
    gs = GridSpec(
        2, 2,
        width_ratios=[1, 5],
        height_ratios=[1, 5],
        hspace=0.02,
        wspace=0.02,
        figure=fig,
    )

    ax_col_dendro = fig.add_subplot(gs[0, 1])
    ax_row_dendro = fig.add_subplot(gs[1, 0])
    ax_heatmap = fig.add_subplot(gs[1, 1])

    # Explicitly pin each dendrogram axis' limits to exactly span its leaves
    # (no autoscale margin). Without this, scipy's leaf spacing (5, 15, 25, ...)
    # and matplotlib's default 5% autoscale margin cause the dendrogram to
    # drift out of alignment with the heatmap columns/rows below/beside it.
    dendrogram(col_linkage, ax=ax_col_dendro, no_labels=True, color_threshold=0)
    ax_col_dendro.set_xlim(0, 10 * n_cols)
    ax_col_dendro.axis("off")

    dendrogram(row_linkage, ax=ax_row_dendro, orientation="left", no_labels=True, color_threshold=0)
    ax_row_dendro.set_ylim(0, 10 * n_rows)
    ax_row_dendro.axis("off")

    im = ax_heatmap.imshow(ordered_matrix, aspect="auto", cmap="YlOrRd")

    ax_heatmap.set_xticks(range(n_cols))
    ax_heatmap.set_xticklabels(ordered_samples, rotation=45, ha="right")
    ax_heatmap.set_yticks(range(n_rows))
    ax_heatmap.set_yticklabels(ordered_species, fontsize=8)
    ax_heatmap.yaxis.tick_right()

    cbar = fig.colorbar(im, ax=ax_heatmap, fraction=0.025, pad=0.35)
    cbar.set_label("log10(CPM + 1)")

    fig.suptitle(f"Top {n_rows} Fungal Species — CPM Heatmap (clustered)", y=0.98)
    fig.savefig(HEATMAP_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Clustered heatmap saved to: {HEATMAP_OUT}")


# ============================================================
# MAIN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    samples = SAMPLES
    if samples is None:
        samples = sorted(
            p.name.replace("_bracken_report.txt", "")
            for p in BRACKEN_REPORTS_DIR.glob("*_bracken_report.txt")
        )
    print(f"Detected samples: {samples}\n")

    print("=" * 60)
    print("STEP 1: Parsing fungal species from each report")
    print("=" * 60)
    all_species = {}
    all_root_totals = {}
    for sample in samples:
        report_path = BRACKEN_REPORTS_DIR / f"{sample}_bracken_report.txt"
        if not report_path.exists():
            print(f"  WARNING: missing {report_path}")
            continue
        root_total, species_counts = parse_fungi_species(report_path)
        all_species[sample] = species_counts
        all_root_totals[sample] = root_total
        print(f"  {sample}: {len(species_counts)} fungal species, root={root_total}")
    print()

    print("=" * 60)
    print(f"STEP 2: Selecting top {TOP_N} species by weighted_abundance_score")
    print("=" * 60)
    top_species = rank_species_by_weighted_score(all_species, all_root_totals, samples, TOP_N)
    for sp in top_species:
        print(f"  {sp}")
    print()

    print("=" * 60)
    print("STEP 3-4: Building matrix and drawing clustered heatmap")
    print("=" * 60)
    matrix = build_matrix(top_species, samples, all_species, all_root_totals)
    plot_clustered_heatmap(matrix, top_species, samples)


if __name__ == "__main__":
    main()
    
    
    