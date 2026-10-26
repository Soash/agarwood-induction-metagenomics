#!/usr/bin/env python3
"""
fungi_summary.py

Fungi-only breakdown: extracts every species-level taxon under kingdom
Fungi (tax_id 4751) from each sample's Bracken report, then produces
CSVs and charts summarizing fungal species composition and abundance.

Input: {sample}_bracken_report.txt files (standard Kraken2/Bracken report
format: percentage, count_clade, count_direct, rank_code, tax_id, name,
with name indented to reflect tree depth).

Method: walks each report top-to-bottom tracking indentation. Once the
Fungi line (tax_id 4751) is reached, every subsequent line with greater
indentation belongs to its subtree; the walk stops at the first line
whose indentation returns to <= the Fungi line's indentation. Within
that subtree, only rank_code == 'S' (species) lines are kept, so counts
match the species-level flat {sample}_bracken.txt tables. Species in
EXCLUDED_SPECIES are dropped entirely at this parsing step, so they
never appear in any downstream CSV, chart, or total.

CPM = raw_species_count / total_classified_reads_in_sample * 1,000,000
(same root-total denominator used in domain_summary.py, tax_id 1).

Outputs (written to OUTPUT_DIR):
    fungi_species_counts.csv               top N species x sample, raw + CPM, "Other fungi" + TOTAL
    fungi_species_presence_cpm.csv         ALL species: presence/prevalence/CPM/weighted_abundance_score
    fungi_species_cpm_by_species.png       grouped bar, x=species (top N), legend=sample (log scale)
    fungi_species_cpm_by_sample.png        grouped bar, x=sample, legend=species (top N) (log scale)
    fungi_species_abundance_by_sample.png  grouped bar, x=sample, legend=species (top N by weighted score)
    fungi_species_richness_by_sample.png   number of distinct fungal species detected per sample
    fungi_composition_stackedbar.png       100%-stacked bar, fungal species composition per sample
    fungi_total_cpm_barchart.png           total fungal CPM per sample (single bar per sample)

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/fungi_summary.py
"""

import csv
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
except ImportError:
    plt = None

# ============================================================
# CONFIG — adjust these if your layout differs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/

BRACKEN_REPORTS_DIR = PROJECT_ROOT / "bracken_results"
OUTPUT_DIR = PROJECT_ROOT / "fungi_summary"

CSV_OUT = OUTPUT_DIR / "fungi_species_counts.csv"
PRESENCE_CPM_CSV_OUT = OUTPUT_DIR / "fungi_species_presence_cpm.csv"
GROUPEDBAR_BYSPECIES_OUT = OUTPUT_DIR / "fungi_species_cpm_by_species.png"
GROUPEDBAR_BYSAMPLE_OUT = OUTPUT_DIR / "fungi_species_cpm_by_sample.png"
STACKEDBAR_OUT = OUTPUT_DIR / "fungi_composition_stackedbar.png"
TOTAL_CPM_BAR_OUT = OUTPUT_DIR / "fungi_total_cpm_barchart.png"
ABUNDANCE_BYSAMPLE_OUT = OUTPUT_DIR / "fungi_species_abundance_by_sample.png"
RICHNESS_BARCHART_OUT = OUTPUT_DIR / "fungi_species_richness_by_sample.png"

TAXID_FUNGI = 4751
TAXID_ROOT = 1

TOP_N = 10  # number of individual species to show in charts; rest collapse into "Other fungi"

# 0-1, how much weighted_abundance_score favors prevalence (how widespread a
# species is) over raw abundance (how high its CPM is when present).
# 0.5 = equal weight.
PREVALENCE_WEIGHT = 0.5

# Name of the negative-control sample, used for the presence/CPM table below.
# Must match one of the sample names auto-detected from bracken_results/.
CONTROL_SAMPLE = "Control"

EXCLUDED_SPECIES = {                   # species names to drop entirely from the output
    "Malassezia globosa",              # common human skin fungus, likely handler/lab contaminant
    "Malassezia restricta",            # common human skin fungus, likely handler/lab contaminant
    "Malassezia japonica",             # common human skin fungus, likely handler/lab contaminant
    "uncultured fungus",               # unresolved/ambiguous classification
    "uncultured Ascomycota",           # unresolved/ambiguous classification
    "fungal sp.",                      # unresolved/ambiguous classification
}

# Samples are auto-detected from bracken_results/*_bracken_report.txt,
# but you can hardcode instead if you prefer:
# SAMPLES = ["AI_1", "AI_2", "AI_3", "KNI", "NI_1", "NI_2", "Control"]
SAMPLES = None  # None = auto-detect


# ============================================================
# STEP 1: Parse one report -> (root_total, {species_name: count})
# ============================================================

def parse_fungi_species(path: Path):
    """Return (root_total, {species_name: count_clade}) for one report."""
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

            # Already inside the Fungi subtree
            if indent <= fungi_indent:
                within_fungi = False  # exited subtree; keep scanning in case root line comes later
                continue

            if rank_code == "S":
                if name in EXCLUDED_SPECIES:
                    continue
                species_counts[name] = species_counts.get(name, 0) + count_clade

    return root_total, species_counts


# ============================================================
# STEP 2: Pick top N species by total CPM across all samples
# ============================================================

def pick_top_species(all_species, all_root_totals, samples, top_n):
    total_cpm = {}
    for sample in samples:
        root = all_root_totals.get(sample, 0)
        if not root:
            continue
        for species, raw in all_species.get(sample, {}).items():
            total_cpm[species] = total_cpm.get(species, 0) + raw / root * 1_000_000

    ranked = sorted(total_cpm.items(), key=lambda kv: -kv[1])
    top_species = [s for s, _ in ranked[:top_n]]
    other_species = [s for s, _ in ranked[top_n:]]
    return top_species, other_species


# ============================================================
# STEP 2b: Rank species by weighted_abundance_score (prevalence + abundance)
# ============================================================

def rank_species_by_weighted_score(all_species, all_root_totals, samples, top_n):
    """Same scoring logic as write_presence_cpm_csv's weighted_abundance_score,
    used here to pick the top N species for the abundance-by-sample chart."""
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
# STEP 3: Write combined CSV (top N species + "Other fungi" + TOTAL)
# ============================================================

def write_combined_csv(all_species, all_root_totals, samples, top_species, other_species):
    with open(CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Species"]
        for sample in samples:
            header += [f"{sample}_raw", f"{sample}_CPM"]
        writer.writerow(header)

        for species in top_species:
            row = [species]
            for sample in samples:
                raw = all_species.get(sample, {}).get(species, 0)
                root = all_root_totals.get(sample, 0)
                cpm = round(raw / root * 1_000_000, 3) if root else "NA"
                row += [raw, cpm]
            writer.writerow(row)

        if other_species:
            row = [f"Other fungi ({len(other_species)} species)"]
            for sample in samples:
                raw_sum = sum(all_species.get(sample, {}).get(s, 0) for s in other_species)
                root = all_root_totals.get(sample, 0)
                cpm = round(raw_sum / root * 1_000_000, 3) if root else "NA"
                row += [raw_sum, cpm]
            writer.writerow(row)

        totals_row = ["TOTAL fungi"]
        for sample in samples:
            raw_total = sum(all_species.get(sample, {}).values())
            root = all_root_totals.get(sample, 0)
            cpm_total = round(raw_total / root * 1_000_000, 3) if root else "NA"
            totals_row += [raw_total, cpm_total]
        writer.writerow(totals_row)


# ============================================================
# STEP 3b: Write presence/prevalence/CPM table for ALL species
# ============================================================

def write_presence_cpm_csv(all_species, all_root_totals, samples):
    """One row per species (every species seen in any non-control sample),
    with presence across the non-control samples, prevalence count,
    whether it also appears in the control, per-sample CPM (or 'Absent'),
    mean CPM across only the samples where it's present, and a combined
    weighted_abundance_score balancing prevalence vs. abundance."""
    non_control_samples = [s for s in samples if s != CONTROL_SAMPLE]
    if CONTROL_SAMPLE not in samples:
        print(f"  WARNING: control sample '{CONTROL_SAMPLE}' not found among "
              f"detected samples {samples}; present_in_control will be 'No' for everything.")

    all_species_names = set()
    for sample in non_control_samples:
        all_species_names.update(all_species.get(sample, {}).keys())

    rows = []
    for species in all_species_names:
        present_samples = []
        cpm_by_sample = {}
        for sample in non_control_samples:
            raw = all_species.get(sample, {}).get(species, 0)
            root = all_root_totals.get(sample, 0)
            if raw > 0 and root:
                cpm = raw / root * 1_000_000
                cpm_by_sample[sample] = cpm
                present_samples.append(sample)
            else:
                cpm_by_sample[sample] = None  # Absent

        prevalence = len(present_samples)
        present_in_control = "Yes" if all_species.get(CONTROL_SAMPLE, {}).get(species, 0) > 0 else "No"
        mean_cpm_present = (
            round(sum(cpm_by_sample[s] for s in present_samples) / prevalence, 4)
            if prevalence else "NA"
        )

        rows.append({
            "species": species,
            "presence_samples": ";".join(present_samples),
            "prevalence": prevalence,
            "present_in_control": present_in_control,
            "cpm_by_sample": cpm_by_sample,
            "mean_cpm_present_only": mean_cpm_present,
        })

    numeric_means = [
        r["mean_cpm_present_only"] for r in rows
        if isinstance(r["mean_cpm_present_only"], (int, float))
    ]
    max_present_only = max(numeric_means, default=0.0)

    for r in rows:
        prevalence_fraction = r["prevalence"] / len(non_control_samples) if non_control_samples else 0.0
        mean_val = r["mean_cpm_present_only"]
        mean_val = mean_val if isinstance(mean_val, (int, float)) else 0.0
        abundance_fraction = (mean_val / max_present_only) if max_present_only > 0 else 0.0
        score = (PREVALENCE_WEIGHT * prevalence_fraction) + ((1 - PREVALENCE_WEIGHT) * abundance_fraction)
        r["weighted_abundance_score"] = round(score, 4)

    rows.sort(key=lambda r: (-r["weighted_abundance_score"], r["species"]))

    with open(PRESENCE_CPM_CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["species", "presence_samples", "prevalence", "present_in_control"]
        header += [f"cpm_{s}" for s in non_control_samples]
        header += ["mean_cpm_present_only", "weighted_abundance_score"]
        writer.writerow(header)

        for r in rows:
            row = [r["species"], r["presence_samples"], r["prevalence"], r["present_in_control"]]
            for sample in non_control_samples:
                cpm = r["cpm_by_sample"][sample]
                row.append("Absent" if cpm is None else f"{cpm:.4f}")
            mean_val = r["mean_cpm_present_only"]
            row.append(mean_val if mean_val == "NA" else f"{mean_val:.4f}")
            row.append(f"{r['weighted_abundance_score']:.4f}")
            writer.writerow(row)


# ============================================================
# STEP 4a: Grouped bar chart of CPM, x=species (top N), legend=sample
# ============================================================

def plot_cpm_by_species(all_species, all_root_totals, samples, top_species):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping species bar chart.")
        return
    if not top_species:
        print("  WARNING: no fungal species found, skipping species bar chart.")
        return

    n_species = len(top_species)
    n_samples = len(samples)
    x = range(n_species)
    width = 0.8 / n_samples

    fig, ax = plt.subplots(figsize=(12, 6))
    for i, sample in enumerate(samples):
        root = all_root_totals.get(sample, 0)
        values = [
            (all_species.get(sample, {}).get(sp, 0) / root * 1_000_000) if root else 0
            for sp in top_species
        ]
        offsets = [xi + i * width - 0.4 + width / 2 for xi in x]
        ax.bar(offsets, values, width=width, label=sample)

    ax.set_xticks(list(x))
    ax.set_xticklabels(top_species, rotation=45, ha="right")
    ax.set_ylabel("CPM (counts per million classified reads, log scale)")
    ax.set_xlabel("Fungal species")
    ax.set_yscale("log")
    ax.set_title(f"Top {n_species} Fungal Species — CPM by Species and Sample")
    ax.legend(title="Sample", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(GROUPEDBAR_BYSPECIES_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  CPM-by-species bar chart saved to: {GROUPEDBAR_BYSPECIES_OUT}")


# ============================================================
# STEP 4b: Grouped bar chart of CPM, x=sample, legend=species (top N)
# ============================================================

def plot_cpm_by_sample(all_species, all_root_totals, samples, top_species):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping sample bar chart.")
        return
    if not top_species:
        print("  WARNING: no fungal species found, skipping sample bar chart.")
        return

    n_samples = len(samples)
    n_species = len(top_species)
    x = range(n_samples)
    width = 0.8 / n_species

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = cm.tab10.colors if n_species <= 10 else cm.tab20.colors
    for j, species in enumerate(top_species):
        values = []
        for sample in samples:
            root = all_root_totals.get(sample, 0)
            raw = all_species.get(sample, {}).get(species, 0)
            values.append(raw / root * 1_000_000 if root else 0)
        offsets = [xi + j * width - 0.4 + width / 2 for xi in x]
        ax.bar(offsets, values, width=width, label=species, color=colors[j % len(colors)])

    ax.set_xticks(list(x))
    ax.set_xticklabels(samples)
    ax.set_ylabel("CPM (counts per million classified reads, log scale)")
    ax.set_xlabel("Sample")
    ax.set_yscale("log")
    ax.set_title(f"Top {n_species} Fungal Species — CPM by Sample and Species")
    ax.legend(title="Species", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(GROUPEDBAR_BYSAMPLE_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  CPM-by-sample bar chart saved to: {GROUPEDBAR_BYSAMPLE_OUT}")


# ============================================================
# STEP 4b2: Bar graph of species abundance x sample
# (top species ranked by weighted_abundance_score, not just raw CPM)
# ============================================================

def plot_species_abundance_by_sample(all_species, all_root_totals, samples, top_species_weighted):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping species-abundance-by-sample chart.")
        return
    if not top_species_weighted:
        print("  WARNING: no fungal species found, skipping species-abundance-by-sample chart.")
        return

    n_samples = len(samples)
    n_species = len(top_species_weighted)
    x = range(n_samples)
    width = 0.8 / n_species

    fig, ax = plt.subplots(figsize=(11, 6))
    colors = cm.tab10.colors if n_species <= 10 else cm.tab20.colors
    for j, species in enumerate(top_species_weighted):
        values = []
        for sample in samples:
            root = all_root_totals.get(sample, 0)
            raw = all_species.get(sample, {}).get(species, 0)
            values.append(raw / root * 1_000_000 if root else 0)
        offsets = [xi + j * width - 0.4 + width / 2 for xi in x]
        ax.bar(offsets, values, width=width, label=species, color=colors[j % len(colors)])

    ax.set_xticks(list(x))
    ax.set_xticklabels(samples)
    ax.set_ylabel("CPM (counts per million classified reads, log scale)")
    ax.set_xlabel("Sample")
    ax.set_yscale("log")
    ax.set_title(f"Fungal Species Abundance by Sample\n(top {n_species} species by weighted_abundance_score)")
    ax.legend(title="Species", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(ABUNDANCE_BYSAMPLE_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Species-abundance-by-sample bar chart saved to: {ABUNDANCE_BYSAMPLE_OUT}")


# ============================================================
# STEP 4b3: Bar chart of fungal species RICHNESS per sample
# (how many distinct species detected, not their abundance)
# ============================================================

def plot_species_richness_by_sample(all_species, samples):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping species richness chart.")
        return

    values = [len(all_species.get(sample, {})) for sample in samples]

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(samples, values, color=cm.tab10.colors[: len(samples)])
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{val}",
            ha="center",
            va="bottom",
        )

    ax.set_ylabel("Number of distinct fungal species detected")
    ax.set_xlabel("Sample")
    ax.set_title("Fungal Species Richness by Sample")
    fig.savefig(RICHNESS_BARCHART_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Species richness bar chart saved to: {RICHNESS_BARCHART_OUT}")


# ============================================================
# STEP 4c: 100%-stacked bar chart of fungal composition per sample
# ============================================================

def plot_composition_stackedbar(all_species, all_root_totals, samples, top_species, other_species):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping composition stacked bar.")
        return
    if not top_species:
        print("  WARNING: no fungal species found, skipping composition stacked bar.")
        return

    categories = top_species + (["Other fungi"] if other_species else [])
    n_cat = len(categories)
    colors = cm.tab10.colors if n_cat <= 10 else cm.tab20.colors

    fig, ax = plt.subplots(figsize=(8, 6))
    bottoms = [0.0] * len(samples)

    for k, category in enumerate(categories):
        pct_values = []
        for sample in samples:
            fungi_total = sum(all_species.get(sample, {}).values())
            if not fungi_total:
                pct_values.append(0)
                continue
            if category == "Other fungi":
                raw = sum(all_species.get(sample, {}).get(s, 0) for s in other_species)
            else:
                raw = all_species.get(sample, {}).get(category, 0)
            pct_values.append(raw / fungi_total * 100)
        ax.bar(
            samples,
            pct_values,
            bottom=bottoms,
            label=category,
            color=colors[k % len(colors)],
        )
        bottoms = [b + p for b, p in zip(bottoms, pct_values)]

    ax.set_ylabel("Percentage of fungal reads (%)")
    ax.set_xlabel("Sample")
    ax.set_title("Fungal Species Composition by Sample")
    ax.set_ylim(0, 100)
    ax.legend(title="Species", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(STACKEDBAR_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Composition stacked bar chart saved to: {STACKEDBAR_OUT}")


# ============================================================
# STEP 4d: Total fungal CPM per sample (single bar per sample)
# ============================================================

def plot_total_fungi_cpm(all_species, all_root_totals, samples):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping total fungal CPM bar chart.")
        return

    values = []
    for sample in samples:
        raw_total = sum(all_species.get(sample, {}).values())
        root = all_root_totals.get(sample, 0)
        values.append(raw_total / root * 1_000_000 if root else 0)

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(samples, values, color=cm.tab10.colors[: len(samples)])
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{val:.0f}",
            ha="center",
            va="bottom",
        )

    ax.set_ylabel("Total fungal CPM")
    ax.set_xlabel("Sample")
    ax.set_title("Total Fungal CPM by Sample")
    fig.savefig(TOTAL_CPM_BAR_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Total fungal CPM bar chart saved to: {TOTAL_CPM_BAR_OUT}")


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
        fungi_total = sum(species_counts.values())
        print(f"  {sample}: {len(species_counts)} fungal species, "
              f"{fungi_total} fungal reads, root={root_total}")
    print()

    print("=" * 60)
    print(f"STEP 2: Selecting top {TOP_N} fungal species by total CPM")
    print("=" * 60)
    top_species, other_species = pick_top_species(all_species, all_root_totals, samples, TOP_N)
    for sp in top_species:
        print(f"  {sp}")
    if other_species:
        print(f"  ...and {len(other_species)} more collapsed into 'Other fungi'")
    print()

    print("=" * 60)
    print("STEP 3: Writing combined CSV")
    print("=" * 60)
    write_combined_csv(all_species, all_root_totals, samples, top_species, other_species)
    print(f"  CSV written to: {CSV_OUT}")
    print()

    print("=" * 60)
    print("STEP 3b: Writing presence/prevalence/CPM table (all species)")
    print("=" * 60)
    write_presence_cpm_csv(all_species, all_root_totals, samples)
    print(f"  CSV written to: {PRESENCE_CPM_CSV_OUT}")
    print()

    print("=" * 60)
    print("STEP 4: Generating charts")
    print("=" * 60)
    plot_cpm_by_species(all_species, all_root_totals, samples, top_species)
    plot_cpm_by_sample(all_species, all_root_totals, samples, top_species)
    top_species_weighted = rank_species_by_weighted_score(all_species, all_root_totals, samples, TOP_N)
    plot_species_abundance_by_sample(all_species, all_root_totals, samples, top_species_weighted)
    plot_species_richness_by_sample(all_species, samples)
    plot_composition_stackedbar(all_species, all_root_totals, samples, top_species, other_species)
    plot_total_fungi_cpm(all_species, all_root_totals, samples)
    print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for sample in samples:
        fungi_total = sum(all_species.get(sample, {}).values())
        root = all_root_totals.get(sample, "?")
        print(f"  {sample}: {fungi_total} fungal reads / {root} total classified reads")


if __name__ == "__main__":
    main()