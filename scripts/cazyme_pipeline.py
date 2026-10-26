#!/usr/bin/env python3
"""
cazyme_pipeline.py

End-to-end CAZyme class-level abundance pipeline.

Steps:
  1. Load total classified reads per sample from general_stats_with_classification.csv
     (column: classified_reads) -- used as the CPM denominator
  2. Filter each sample's dbCAN/results/{sample}_all_vs_CAZy.tsv by
     identity >= MIN_IDENTITY and alignment length >= MIN_ALNLEN
  3. Extract CAZyme CLASS (GH/GT/PL/CE/AA/CBM/SLH) from the sseqid field
     (first match only, to avoid double-counting repeated tokens in headers)
  4. Tally class counts per sample
  5. Normalize to CPM = (raw / classified_reads) * 1,000,000
  6. Write one merged CSV: Class x (raw + CPM) per sample
  7. Plot a pie chart (with legend + %) and a bar chart of total CPM per sample

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/cazyme_pipeline.py
"""

import csv
import re
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

# ============================================================
# CONFIG — adjust these if your layout differs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/

DBCAN_RESULTS_DIR = PROJECT_ROOT / "dbCAN" / "results"
FILTERED_DIR = DBCAN_RESULTS_DIR / "filtered"
CLASS_COUNTS_DIR = DBCAN_RESULTS_DIR / "class_counts"

# Path to the classification stats CSV (Sample, Raw, Trimmed, Host-Removed, classified_reads)
CLASSIFICATION_STATS_CSV = PROJECT_ROOT / "results" /"general_stats_2.csv"

# kraken2 was run WITH --paired mode, so classified_reads (from Bracken) is a
# PAIR count. Diamond hits are counted per individual read, so we double
# classified_reads here to keep numerator/denominator units consistent.
CLASSIFIED_READS_ARE_PAIRS = True

MERGED_CSV_OUT = PROJECT_ROOT / "dbCAN" / "combined_class_abundance_cpm.csv"
PIECHART_OUT = PROJECT_ROOT / "dbCAN" / "total_cpm_by_sample_piechart.png"
BARCHART_OUT = PROJECT_ROOT / "dbCAN" / "total_cpm_by_sample_barchart.png"

# Filtering thresholds (column 3 = %identity, column 4 = alignment length, per outfmt 6)
MIN_IDENTITY = 60.0
MIN_ALNLEN = 30

# CAZyme class regex, in preferred display order
CLASS_ORDER = ["GH", "GT", "PL", "CE", "AA", "CBM", "SLH"]
CLASS_REGEX = re.compile(r"(GH|GT|PL|CE|CBM|AA|SLH)")

# Samples are auto-detected from dbCAN/results/*_all_vs_CAZy.tsv,
# but you can hardcode instead if you prefer:
# SAMPLES = ["AI_1", "AI_2", "AI_3", "KNI", "NI_1", "NI_2", "Control"]
SAMPLES = None  # None = auto-detect


# ============================================================
# STEP 1: Load classified read counts (CPM denominator)
# ============================================================

def load_classified_reads(csv_path: Path):
    """Read Sample -> classified_reads from the stats CSV."""
    totals = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["Sample"].strip()
            try:
                n = int(row["classified_reads"])
                if CLASSIFIED_READS_ARE_PAIRS:
                    n *= 2
                totals[sample] = n
            except (KeyError, ValueError):
                print(f"  WARNING: could not parse classified_reads for row: {row}")
    return totals


# ============================================================
# STEP 2: Filter diamond hits by identity + alignment length
# ============================================================

def filter_sample(sample):
    """Filter {sample}_all_vs_CAZy.tsv, write filtered file, return line count."""
    in_path = DBCAN_RESULTS_DIR / f"{sample}_all_vs_CAZy.tsv"
    out_path = FILTERED_DIR / f"{sample}_all_vs_CAZy.filtered.tsv"

    if not in_path.exists():
        print(f"  WARNING: missing {in_path}")
        return 0

    kept = 0
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                identity = float(fields[2])
                alnlen = float(fields[3])
            except ValueError:
                continue
            if identity >= MIN_IDENTITY and alnlen >= MIN_ALNLEN:
                fout.write(line)
                kept += 1

    return kept


# ============================================================
# STEP 3 + 4: Extract CAZyme class per hit, tally per sample
# ============================================================

def extract_first_class(sseqid: str):
    """Return the first CAZyme class token found in sseqid, or None."""
    m = CLASS_REGEX.search(sseqid)
    return m.group(1) if m else None


def class_count_sample(sample):
    """Read filtered tsv, extract class per line, tally, write count file."""
    in_path = FILTERED_DIR / f"{sample}_all_vs_CAZy.filtered.tsv"
    out_path = CLASS_COUNTS_DIR / f"{sample}_class_counts.txt"

    counts = {}

    if not in_path.exists():
        print(f"  WARNING: missing {in_path}")
        return counts

    with open(in_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            sseqid = fields[1]
            cls = extract_first_class(sseqid)
            if cls:
                counts[cls] = counts.get(cls, 0) + 1

    # Write human-readable count file (same style as `uniq -c | sort -rn`)
    with open(out_path, "w") as f:
        for cls, n in sorted(counts.items(), key=lambda x: -x[1]):
            f.write(f"{n:>7} {cls}\n")

    return counts


# ============================================================
# STEP 5 + 6: Compute CPM and write merged CSV
# ============================================================

def write_merged_csv(all_counts, classified_reads, samples):
    seen_classes = set()
    for counts in all_counts.values():
        seen_classes.update(counts.keys())

    classes = [c for c in CLASS_ORDER if c in seen_classes]
    classes += sorted(seen_classes - set(CLASS_ORDER))

    header = ["Class"]
    for sample in samples:
        header += [f"{sample}_raw", f"{sample}_CPM"]

    with open(MERGED_CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for cls in classes:
            row = [cls]
            for sample in samples:
                raw = all_counts.get(sample, {}).get(cls, 0)
                depth = classified_reads.get(sample)
                cpm = round(raw / depth * 1_000_000, 2) if depth else "NA"
                row += [raw, cpm]
            writer.writerow(row)

        # Totals row
        totals_row = ["TOTAL"]
        for sample in samples:
            raw_total = sum(all_counts.get(sample, {}).values())
            depth = classified_reads.get(sample)
            cpm_total = round(raw_total / depth * 1_000_000, 2) if depth else "NA"
            totals_row += [raw_total, cpm_total]
        writer.writerow(totals_row)

    return classes


def compute_total_cpm_per_sample(all_counts, classified_reads, samples):
    """Return {sample: total_CPM_across_all_classes}."""
    totals = {}
    for sample in samples:
        raw_total = sum(all_counts.get(sample, {}).values())
        depth = classified_reads.get(sample)
        totals[sample] = round(raw_total / depth * 1_000_000, 2) if depth else 0
    return totals


def plot_total_cpm_piechart(total_cpm_per_sample, samples):
    """Pie chart of each sample's share of total CPM (summed across all CAZyme classes).
    No text inside the wedges; a legend on the side shows Sample (percentage%)."""
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping pie chart. "
              "Install with: pip install matplotlib")
        return

    values = [total_cpm_per_sample.get(s, 0) for s in samples]
    total = sum(values)
    if total == 0:
        print("  WARNING: all total CPM values are zero, skipping pie chart.")
        return

    fig, ax = plt.subplots(figsize=(7, 5.5))
    wedges, _ = ax.pie(values, startangle=90)

    legend_labels = [
        f"{s} ({v/total*100:.1f}%)" for s, v in zip(samples, values)
    ]
    ax.legend(
        wedges,
        legend_labels,
        title="Sample",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
    )
    ax.set_title("Total CAZyme CPM by Sample")
    ax.axis("equal")
    fig.savefig(PIECHART_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Pie chart saved to: {PIECHART_OUT}")

def plot_total_cpm_barchart(total_cpm_per_sample, samples):
    """Bar chart of total CPM (summed across all CAZyme classes) per sample."""
    if plt is None:
        print("  WARNING: matplotlib not installed, skipping bar chart. "
              "Install with: pip install matplotlib")
        return

    values = [total_cpm_per_sample.get(s, 0) for s in samples]
    if sum(values) == 0:
        print("  WARNING: all total CPM values are zero, skipping bar chart.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(samples, values, color=plt.cm.tab10.colors[:len(samples)])

    # Label each bar with its CPM value on top
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{val:.0f}",
            ha="center",
            va="bottom",
        )

    ax.set_xlabel("Sample")
    ax.set_ylabel("Total CPM")
    ax.set_title("Total CAZyme CPM by Sample")
    fig.tight_layout()
    fig.savefig(BARCHART_OUT, dpi=200)
    plt.close(fig)
    print(f"  Bar chart saved to: {BARCHART_OUT}")


# ============================================================
# MAIN
# ============================================================

def main():
    FILTERED_DIR.mkdir(parents=True, exist_ok=True)
    CLASS_COUNTS_DIR.mkdir(parents=True, exist_ok=True)

    samples = SAMPLES
    if samples is None:
        samples = sorted(
            p.name.replace("_all_vs_CAZy.tsv", "")
            for p in DBCAN_RESULTS_DIR.glob("*_all_vs_CAZy.tsv")
        )
    print(f"Detected samples: {samples}\n")

    print("=" * 60)
    print("STEP 1: Loading classified_reads (CPM denominator)")
    print("=" * 60)
    if not CLASSIFICATION_STATS_CSV.exists():
        raise FileNotFoundError(
            f"Could not find {CLASSIFICATION_STATS_CSV}. "
            f"Update CLASSIFICATION_STATS_CSV at the top of this script."
        )
    classified_reads = load_classified_reads(CLASSIFICATION_STATS_CSV)
    for sample in samples:
        val = classified_reads.get(sample, "MISSING")
        print(f"  {sample}: classified_reads={val}")
    print()

    print("=" * 60)
    print(f"STEP 2: Filtering hits (identity >= {MIN_IDENTITY}, alnlen >= {MIN_ALNLEN})")
    print("=" * 60)
    filtered_counts = {}
    for sample in samples:
        n = filter_sample(sample)
        filtered_counts[sample] = n
        print(f"  {sample}: {n} hits kept")
    print()

    print("=" * 60)
    print("STEP 3-4: Extracting CAZyme classes and tallying per sample")
    print("=" * 60)
    all_counts = {}
    for sample in samples:
        counts = class_count_sample(sample)
        all_counts[sample] = counts
        classified_total = sum(counts.values())
        filt_total = filtered_counts.get(sample, 0)
        flag = "" if classified_total == filt_total else "  <-- MISMATCH, check regex/headers"
        print(f"  {sample}: {classified_total} classified / {filt_total} filtered hits{flag}")
    print()

    print("=" * 60)
    print("STEP 5-6: Computing CPM and writing merged CSV")
    print("=" * 60)
    write_merged_csv(all_counts, classified_reads, samples)
    print(f"  Merged CSV written to: {MERGED_CSV_OUT}")
    print()

    print("=" * 60)
    print("STEP 7: Plotting total CPM per sample (pie + bar charts)")
    print("=" * 60)
    total_cpm_per_sample = compute_total_cpm_per_sample(all_counts, classified_reads, samples)
    for sample in samples:
        print(f"  {sample}: total CPM = {total_cpm_per_sample.get(sample, 0)}")
    plot_total_cpm_piechart(total_cpm_per_sample, samples)
    plot_total_cpm_barchart(total_cpm_per_sample, samples)
    print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for sample in samples:
        depth = classified_reads.get(sample, "?")
        n_classified = sum(all_counts.get(sample, {}).values())
        print(f"  {sample}: {n_classified} classified CAZyme hits / {depth} classified_reads")


if __name__ == "__main__":
    main()