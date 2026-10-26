#!/usr/bin/env python3
"""
fungi_AI_vs_NI_stats.py

Tests whether fungal species differ significantly in PRESENCE or ABUNDANCE
between two sample groups:
    Artificially Induced (AI): AI_1, AI_2, AI_3
    Naturally Induced   (NI): KNI, NI_1, NI_2

For each fungal species detected in either group (Control is excluded from
this comparison entirely -- it's a negative control, not part of either
biological group), this script runs two independent tests:

  1. Fisher's exact test on presence/absence
     2x2 table: [present_in_AI, absent_in_AI] x [present_in_NI, absent_in_NI]
     Tests whether the species' presence/absence pattern is associated
     with group membership.

  2. Mann-Whitney U test on abundance (CPM)
     Compares the 3 AI CPM values against the 3 NI CPM values (absent
     samples contribute CPM=0, since absence is itself informative for
     an abundance comparison). Non-parametric, doesn't assume normality.

Both p-values are FDR-corrected (Benjamini-Hochberg) across all species
tested, since many species are tested simultaneously.

*** IMPORTANT CAVEAT ***
With only 3 samples per group, statistical power is very low. Fisher's
exact test with n=3 vs n=3 can only take a handful of possible p-values
(the smallest achievable two-sided p-value is 0.1, and only for a perfect
3-vs-0 split), and Mann-Whitney U similarly has very limited resolution.
Treat any results here as exploratory / hypothesis-generating, not
confirmatory. Do not over-interpret non-significant results as "no
difference" -- absence of significance with n=3 per group is expected
even for real, meaningful differences.

Uses the same Fungi-subtree parsing and EXCLUDED_SPECIES contaminant list
as fungi_summary.py, so results are consistent between the two scripts.

Outputs (written to OUTPUT_DIR):
    fungi_AI_vs_NI_stats.csv         full results table, all species tested
    fungi_AI_vs_NI_volcano.png       log2FC (mean CPM) vs -log10(Mann-Whitney p)
    fungi_AI_vs_NI_significant_barchart.png   count of species enriched in each group

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/fungi_AI_vs_NI_stats.py
"""

import csv
from pathlib import Path

from scipy.stats import fisher_exact, mannwhitneyu

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

# ============================================================
# CONFIG — adjust these if your layout differs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/

BRACKEN_REPORTS_DIR = PROJECT_ROOT / "bracken_results"
OUTPUT_DIR = PROJECT_ROOT / "fungi_summary"

CSV_OUT = OUTPUT_DIR / "fungi_AI_vs_NI_stats.csv"
VOLCANO_OUT = OUTPUT_DIR / "fungi_AI_vs_NI_volcano.png"
SIG_BARCHART_OUT = OUTPUT_DIR / "fungi_AI_vs_NI_significant_barchart.png"

TAXID_FUNGI = 4751
TAXID_ROOT = 1

AI_GROUP = ["AI_1", "AI_2", "AI_3"]     # Artificially Induced
NI_GROUP = ["KNI", "NI_1", "NI_2"]      # Naturally Induced

FDR_ALPHA = 0.05     # significance threshold applied to BH-adjusted q-values
LOG2FC_PSEUDOCOUNT = 1.0  # added to mean CPM before log2 fold-change, avoids log(0)

EXCLUDED_SPECIES = {                   # same contaminant/ambiguous list as fungi_summary.py
    "Malassezia globosa",
    "Malassezia restricta",
    "Malassezia japonica",
    "uncultured fungus",
    "uncultured Ascomycota",
    "fungal sp.",
}


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
# STEP 2: Benjamini-Hochberg FDR correction
# ============================================================

def bh_correction(pvalues):
    """Standard Benjamini-Hochberg step-up procedure. Returns q-values in
    the same order as the input list."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        p = pvalues[idx]
        val = min(p * n / rank, 1.0)
        running_min = min(running_min, val)
        q[idx] = running_min
    return q


# ============================================================
# STEP 3: Run the AI vs NI comparison for one species
# ============================================================

def compare_species(species, all_species, all_root_totals):
    cpm_ai = []
    present_ai = 0
    for sample in AI_GROUP:
        raw = all_species.get(sample, {}).get(species, 0)
        root = all_root_totals.get(sample, 0)
        cpm = (raw / root * 1_000_000) if root else 0.0
        cpm_ai.append(cpm)
        if raw > 0:
            present_ai += 1

    cpm_ni = []
    present_ni = 0
    for sample in NI_GROUP:
        raw = all_species.get(sample, {}).get(species, 0)
        root = all_root_totals.get(sample, 0)
        cpm = (raw / root * 1_000_000) if root else 0.0
        cpm_ni.append(cpm)
        if raw > 0:
            present_ni += 1

    # Fisher's exact test on presence/absence
    table = [
        [present_ai, len(AI_GROUP) - present_ai],
        [present_ni, len(NI_GROUP) - present_ni],
    ]
    try:
        _, fisher_p = fisher_exact(table)
    except ValueError:
        fisher_p = 1.0

    # Mann-Whitney U test on abundance (CPM), zeros included
    try:
        _, mw_p = mannwhitneyu(cpm_ai, cpm_ni, alternative="two-sided")
    except ValueError:
        # raised when all values in both groups are identical (e.g. all zero)
        mw_p = 1.0

    mean_ai = sum(cpm_ai) / len(cpm_ai)
    mean_ni = sum(cpm_ni) / len(cpm_ni)
    log2fc = __import__("math").log2(
        (mean_ai + LOG2FC_PSEUDOCOUNT) / (mean_ni + LOG2FC_PSEUDOCOUNT)
    )

    return {
        "species": species,
        "presence_AI": f"{present_ai}/{len(AI_GROUP)}",
        "presence_NI": f"{present_ni}/{len(NI_GROUP)}",
        "cpm_ai_by_sample": dict(zip(AI_GROUP, cpm_ai)),
        "cpm_ni_by_sample": dict(zip(NI_GROUP, cpm_ni)),
        "mean_cpm_AI": mean_ai,
        "mean_cpm_NI": mean_ni,
        "log2fc_AI_vs_NI": log2fc,
        "fisher_pvalue": fisher_p,
        "mannwhitney_pvalue": mw_p,
    }


# ============================================================
# STEP 4: Write results CSV
# ============================================================

def write_results_csv(results):
    fisher_p = [r["fisher_pvalue"] for r in results]
    mw_p = [r["mannwhitney_pvalue"] for r in results]
    fisher_q = bh_correction(fisher_p)
    mw_q = bh_correction(mw_p)

    for r, fq, mq in zip(results, fisher_q, mw_q):
        r["fisher_qvalue"] = fq
        r["mannwhitney_qvalue"] = mq
        r["significant_presence"] = "Yes" if fq < FDR_ALPHA else "No"
        r["significant_abundance"] = "Yes" if mq < FDR_ALPHA else "No"

    results.sort(key=lambda r: min(r["fisher_qvalue"], r["mannwhitney_qvalue"]))

    with open(CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        header = (
            ["species", "presence_AI", "presence_NI"]
            + [f"cpm_{s}" for s in AI_GROUP]
            + [f"cpm_{s}" for s in NI_GROUP]
            + ["mean_cpm_AI", "mean_cpm_NI", "log2fc_AI_vs_NI",
               "fisher_pvalue", "fisher_qvalue", "significant_presence",
               "mannwhitney_pvalue", "mannwhitney_qvalue", "significant_abundance"]
        )
        writer.writerow(header)

        for r in results:
            row = [r["species"], r["presence_AI"], r["presence_NI"]]
            row += [f"{r['cpm_ai_by_sample'][s]:.4f}" for s in AI_GROUP]
            row += [f"{r['cpm_ni_by_sample'][s]:.4f}" for s in NI_GROUP]
            row += [
                f"{r['mean_cpm_AI']:.4f}",
                f"{r['mean_cpm_NI']:.4f}",
                f"{r['log2fc_AI_vs_NI']:.4f}",
                f"{r['fisher_pvalue']:.4f}",
                f"{r['fisher_qvalue']:.4f}",
                r["significant_presence"],
                f"{r['mannwhitney_pvalue']:.4f}",
                f"{r['mannwhitney_qvalue']:.4f}",
                r["significant_abundance"],
            ]
            writer.writerow(row)


# ============================================================
# STEP 5a: Volcano plot (log2FC vs -log10 Mann-Whitney p)
# ============================================================

def plot_volcano(results):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping volcano plot.")
        return

    import math

    x = [r["log2fc_AI_vs_NI"] for r in results]
    y = [-math.log10(r["mannwhitney_pvalue"]) if r["mannwhitney_pvalue"] > 0 else 0 for r in results]
    colors = ["#C44E52" if r["significant_abundance"] == "Yes" else "#8C8C8C" for r in results]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, c=colors, alpha=0.7, edgecolors="none")
    ax.axhline(-math.log10(0.05), color="black", linestyle="--", linewidth=0.8, label="raw p = 0.05")
    ax.axvline(0, color="black", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log2 fold-change (AI vs NI, mean CPM)")
    ax.set_ylabel("-log10(Mann-Whitney p-value)")
    ax.set_title("Fungal Species: AI vs NI\n(red = significant at FDR < %.2f)" % FDR_ALPHA)
    ax.legend()
    fig.savefig(VOLCANO_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Volcano plot saved to: {VOLCANO_OUT}")


# ============================================================
# STEP 5b: Bar chart summarizing significant hits per direction
# ============================================================

def plot_significant_summary(results):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping significance summary chart.")
        return

    enriched_ai = sum(
        1 for r in results if r["significant_abundance"] == "Yes" and r["log2fc_AI_vs_NI"] > 0
    )
    enriched_ni = sum(
        1 for r in results if r["significant_abundance"] == "Yes" and r["log2fc_AI_vs_NI"] < 0
    )
    presence_sig = sum(1 for r in results if r["significant_presence"] == "Yes")

    categories = ["Enriched in AI\n(abundance)", "Enriched in NI\n(abundance)", "Presence/absence\nsignificant"]
    values = [enriched_ai, enriched_ni, presence_sig]

    fig, ax = plt.subplots(figsize=(7, 6))
    bars = ax.bar(categories, values, color=["#4C72B0", "#DD8452", "#55A868"])
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(val),
                ha="center", va="bottom")

    ax.set_ylabel(f"Number of species (FDR < {FDR_ALPHA})")
    ax.set_title("Significant Fungal Species: AI vs NI")
    fig.savefig(SIG_BARCHART_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Significance summary bar chart saved to: {SIG_BARCHART_OUT}")


# ============================================================
# MAIN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_samples = AI_GROUP + NI_GROUP
    print(f"AI group: {AI_GROUP}")
    print(f"NI group: {NI_GROUP}\n")

    print("=" * 60)
    print("STEP 1: Parsing fungal species from each report")
    print("=" * 60)
    all_species = {}
    all_root_totals = {}
    for sample in all_samples:
        report_path = BRACKEN_REPORTS_DIR / f"{sample}_bracken_report.txt"
        if not report_path.exists():
            print(f"  WARNING: missing {report_path}")
            continue
        root_total, species_counts = parse_fungi_species(report_path)
        all_species[sample] = species_counts
        all_root_totals[sample] = root_total
        print(f"  {sample}: {len(species_counts)} fungal species, root={root_total}")
    print()

    all_species_names = set()
    for sample in all_samples:
        all_species_names.update(all_species.get(sample, {}).keys())
    print(f"Total distinct fungal species across AI+NI: {len(all_species_names)}\n")

    print("=" * 60)
    print("STEP 2-3: Running Fisher's exact + Mann-Whitney U per species")
    print("=" * 60)
    results = [compare_species(sp, all_species, all_root_totals) for sp in all_species_names]
    print(f"  Tested {len(results)} species.\n")

    print("=" * 60)
    print("STEP 4: Writing results CSV")
    print("=" * 60)
    write_results_csv(results)
    print(f"  CSV written to: {CSV_OUT}")
    print()

    print("=" * 60)
    print("STEP 5: Generating charts")
    print("=" * 60)
    plot_volcano(results)
    plot_significant_summary(results)
    print()

    n_sig_presence = sum(1 for r in results if r["significant_presence"] == "Yes")
    n_sig_abundance = sum(1 for r in results if r["significant_abundance"] == "Yes")
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Species with significant presence/absence difference (FDR<{FDR_ALPHA}): {n_sig_presence}")
    print(f"  Species with significant abundance difference (FDR<{FDR_ALPHA}): {n_sig_abundance}")
    print()
    print("  *** CAVEAT: n=3 per group gives very low statistical power. ***")
    print("  *** Treat results as exploratory / hypothesis-generating.  ***")


if __name__ == "__main__":
    main()