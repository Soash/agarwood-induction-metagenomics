#!/usr/bin/env python3
"""
domain_summary.py

Summarize Bracken-corrected read counts by top-level taxonomic domain:
Bacteria, Fungi, Archaea, Virus, and Others (all non-fungal Eukaryota:
plants, animals, protists, etc.)

Input: {sample}_bracken_report.txt files (standard Kraken2/Bracken report
format: percentage, count_clade, count_direct, rank_code, tax_id, name).
These are the tree-structured reports (NOT the flat species-level
{sample}_bracken.txt files), since only the tree report lets us read off
domain-level clade totals directly by tax_id.

Classification uses fixed NCBI taxonomy IDs (stable, rank-independent):
    Bacteria   tax_id 2
    Archaea    tax_id 2157
    Eukaryota  tax_id 2759
    Viruses    tax_id 10239
    Fungi      tax_id 4751   (kingdom, nested inside Eukaryota)
    Others   = Eukaryota_total - Fungi_total  (plants, animals, protists...)

Species in EXCLUDED_SPECIES (e.g. likely handler/lab contaminants) are
subtracted from the Fungi total specifically -- they're dropped entirely,
not reassigned to "Others", since biologically they are fungi, just
excluded as noise. The CPM denominator (root total) is left unchanged,
matching the same convention used in fungi_summary.py.

Outputs (written to OUTPUT_DIR):
    combined_domain_counts.csv           raw + percentage per sample per domain
    domain_composition_stackedbar.png    100%-stacked bar, composition per sample
    domain_cpm_barchart_by_domain.png    grouped bar chart, CPM (log scale), x=Domain, legend=Sample
    domain_cpm_barchart_by_sample.png    grouped bar chart, CPM (log scale), x=Sample, legend=Domain

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/domain_summary.py
"""

import csv
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

# ============================================================
# CONFIG — adjust these if your layout differs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/

BRACKEN_REPORTS_DIR = PROJECT_ROOT / "bracken_results"
OUTPUT_DIR = PROJECT_ROOT / "domain_summary"

CSV_OUT = OUTPUT_DIR / "combined_domain_counts.csv"
STACKEDBAR_OUT = OUTPUT_DIR / "domain_composition_stackedbar.png"
GROUPEDBAR_CPM_BYDOMAIN_OUT = OUTPUT_DIR / "domain_cpm_barchart_by_domain.png"
GROUPEDBAR_CPM_BYSAMPLE_OUT = OUTPUT_DIR / "domain_cpm_barchart_by_sample.png"

# Fixed NCBI taxonomy IDs — stable regardless of rank-code position in the tree
TAXID_BACTERIA = 2
TAXID_ARCHAEA = 2157
TAXID_EUKARYOTA = 2759
TAXID_VIRUSES = 10239
TAXID_FUNGI = 4751
TAXID_ROOT = 1

DOMAIN_ORDER = ["Bacteria", "Fungi", "Archaea", "Virus", "Others"]
DOMAIN_COLORS = {
    "Bacteria": "#4C72B0",
    "Fungi": "#DD8452",
    "Archaea": "#55A868",
    "Virus": "#C44E52",
    "Others": "#8172B2",
}

# Species names to drop from the Fungi total entirely (likely handler/lab
# contaminants or unresolved/ambiguous classifications). Must match the
# 'name' field exactly as it appears in the report (whitespace-stripped).
EXCLUDED_SPECIES = {
    "Malassezia globosa",              # common human skin fungus, likely handler/lab contaminant
    "Malassezia restricta",            # common human skin fungus, likely handler/lab contaminant
    "Malassezia japonica",             # common human skin fungus, likely handler/lab contaminant
    "uncultured fungus",               # unresolved/ambiguous classification
    "uncultured Ascomycota",           # unresolved/ambiguous classification
    "fungal sp.",           # unresolved/ambiguous classification
}

# Samples are auto-detected from bracken_results/*_bracken_report.txt,
# but you can hardcode instead if you prefer:
# SAMPLES = ["AI_1", "AI_2", "AI_3", "KNI", "NI_1", "NI_2", "Control"]
SAMPLES = None  # None = auto-detect


# ============================================================
# STEP 1: Parse a bracken_report.txt -> clade_counts + excluded total
# ============================================================

def parse_report_clade_counts(path: Path):
    """Return ({tax_id: count_clade}, excluded_species_total) from a
    Kraken2/Bracken-style report. excluded_species_total is the sum of
    count_clade for every species-level (rank 'S') line whose name is
    in EXCLUDED_SPECIES."""
    clade_counts = {}
    excluded_total = 0

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

            clade_counts[tax_id] = count_clade

            if rank_code == "S":
                name = fields[5].strip()
                if name in EXCLUDED_SPECIES:
                    excluded_total += count_clade

    return clade_counts, excluded_total


# ============================================================
# STEP 2: Compute domain-level totals for one sample
# ============================================================

def compute_domain_totals(clade_counts, excluded_total):
    bacteria = clade_counts.get(TAXID_BACTERIA, 0)
    archaea = clade_counts.get(TAXID_ARCHAEA, 0)
    eukaryota = clade_counts.get(TAXID_EUKARYOTA, 0)
    viruses = clade_counts.get(TAXID_VIRUSES, 0)
    fungi = max(clade_counts.get(TAXID_FUNGI, 0) - excluded_total, 0)
    others = max(eukaryota - clade_counts.get(TAXID_FUNGI, 0), 0)  # non-fungal eukaryotes

    totals = {
        "Bacteria": bacteria,
        "Fungi": fungi,
        "Archaea": archaea,
        "Virus": viruses,
        "Others": others,
    }

    root_total = clade_counts.get(TAXID_ROOT, sum(totals.values()))
    return totals, root_total


# ============================================================
# STEP 3: Write combined CSV (raw + percentage per sample)
# ============================================================

def write_combined_csv(all_totals, all_root_totals, samples):
    with open(CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Domain"]
        for sample in samples:
            header += [f"{sample}_raw", f"{sample}_pct"]
        writer.writerow(header)

        for domain in DOMAIN_ORDER:
            row = [domain]
            for sample in samples:
                raw = all_totals.get(sample, {}).get(domain, 0)
                root = all_root_totals.get(sample, 0)
                pct = round(raw / root * 100, 3) if root else "NA"
                row += [raw, pct]
            writer.writerow(row)

        totals_row = ["TOTAL (classified)"]
        for sample in samples:
            root = all_root_totals.get(sample, 0)
            totals_row += [root, 100.0 if root else "NA"]
        writer.writerow(totals_row)


# ============================================================
# STEP 4a: 100%-stacked bar chart (composition per sample)
# ============================================================

def plot_stacked_composition(all_totals, all_root_totals, samples):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping stacked bar chart.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    bottoms = [0.0] * len(samples)

    for domain in DOMAIN_ORDER:
        pct_values = []
        for sample in samples:
            raw = all_totals.get(sample, {}).get(domain, 0)
            root = all_root_totals.get(sample, 0)
            pct_values.append(raw / root * 100 if root else 0)
        ax.bar(
            samples,
            pct_values,
            bottom=bottoms,
            label=domain,
            color=DOMAIN_COLORS.get(domain),
        )
        bottoms = [b + p for b, p in zip(bottoms, pct_values)]

    ax.set_ylabel("Percentage of classified reads (%)")
    ax.set_xlabel("Sample")
    ax.set_title("Taxonomic Domain Composition by Sample")
    ax.set_ylim(0, 100)
    ax.legend(title="Domain", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(STACKEDBAR_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Stacked composition bar chart saved to: {STACKEDBAR_OUT}")


# ============================================================
# STEP 4b: Grouped bar chart of CPM (log scale), x=Domain, legend=Sample
# ============================================================

def plot_grouped_cpm_by_domain(all_totals, all_root_totals, samples):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping CPM-by-domain bar chart.")
        return

    n_domains = len(DOMAIN_ORDER)
    n_samples = len(samples)
    x = range(n_domains)
    width = 0.8 / n_samples

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, sample in enumerate(samples):
        root = all_root_totals.get(sample, 0)
        values = [
            (all_totals.get(sample, {}).get(d, 0) / root * 1_000_000) if root else 0
            for d in DOMAIN_ORDER
        ]
        offsets = [xi + i * width - 0.4 + width / 2 for xi in x]
        ax.bar(offsets, values, width=width, label=sample)

    ax.set_xticks(list(x))
    ax.set_xticklabels(DOMAIN_ORDER)
    ax.set_ylabel("CPM (counts per million classified reads, log scale)")
    ax.set_xlabel("Domain")
    ax.set_yscale("log")
    ax.set_title("CPM by Domain and Sample")
    ax.legend(title="Sample", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(GROUPEDBAR_CPM_BYDOMAIN_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Grouped CPM-by-domain bar chart saved to: {GROUPEDBAR_CPM_BYDOMAIN_OUT}")


# ============================================================
# STEP 4c: Grouped bar chart of CPM, x-axis=Sample, legend=Domain
# ============================================================

def plot_grouped_cpm_by_sample(all_totals, all_root_totals, samples):
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping CPM-by-sample bar chart.")
        return

    n_samples = len(samples)
    n_domains = len(DOMAIN_ORDER)
    x = range(n_samples)
    width = 0.8 / n_domains

    fig, ax = plt.subplots(figsize=(10, 6))
    for j, domain in enumerate(DOMAIN_ORDER):
        values = []
        for sample in samples:
            root = all_root_totals.get(sample, 0)
            raw = all_totals.get(sample, {}).get(domain, 0)
            values.append(raw / root * 1_000_000 if root else 0)
        offsets = [xi + j * width - 0.4 + width / 2 for xi in x]
        ax.bar(offsets, values, width=width, label=domain, color=DOMAIN_COLORS.get(domain))

    ax.set_xticks(list(x))
    ax.set_xticklabels(samples)
    ax.set_ylabel("CPM (counts per million classified reads, log scale)")
    ax.set_xlabel("Sample")
    ax.set_yscale("log")
    ax.set_title("CPM by Sample and Domain")
    ax.legend(title="Domain", loc="center left", bbox_to_anchor=(1, 0.5))
    fig.savefig(GROUPEDBAR_CPM_BYSAMPLE_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Grouped CPM-by-sample bar chart saved to: {GROUPEDBAR_CPM_BYSAMPLE_OUT}")


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
    print("STEP 1-2: Parsing reports and computing domain totals")
    print("=" * 60)
    all_totals = {}
    all_root_totals = {}
    for sample in samples:
        report_path = BRACKEN_REPORTS_DIR / f"{sample}_bracken_report.txt"
        if not report_path.exists():
            print(f"  WARNING: missing {report_path}")
            continue
        clade_counts, excluded_total = parse_report_clade_counts(report_path)
        totals, root_total = compute_domain_totals(clade_counts, excluded_total)
        all_totals[sample] = totals
        all_root_totals[sample] = root_total

        domain_sum = sum(totals.values())
        flag = "" if abs(domain_sum - root_total) <= 1 else "  <-- check: domain sum != root total"
        excl_note = f"  (excluded {excluded_total} reads)" if excluded_total else ""
        print(f"  {sample}: " + ", ".join(f"{d}={totals[d]}" for d in DOMAIN_ORDER) +
              f"  | root={root_total}{flag}{excl_note}")
    print()

    print("=" * 60)
    print("STEP 3: Writing combined CSV")
    print("=" * 60)
    write_combined_csv(all_totals, all_root_totals, samples)
    print(f"  CSV written to: {CSV_OUT}")
    print()

    print("=" * 60)
    print("STEP 4: Generating charts")
    print("=" * 60)
    plot_stacked_composition(all_totals, all_root_totals, samples)
    plot_grouped_cpm_by_domain(all_totals, all_root_totals, samples)
    plot_grouped_cpm_by_sample(all_totals, all_root_totals, samples)
    print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for sample in samples:
        root = all_root_totals.get(sample, "?")
        print(f"  {sample}: {root} total classified reads")


if __name__ == "__main__":
    main()