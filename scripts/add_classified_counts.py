"""add_classified_counts.py

Takes an existing stats CSV (one row per sample, e.g. Raw/Trimmed/Host-Removed
counts from MultiQC) and adds a new column pulled from each sample's
Kraken/Bracken report file:

    classified_reads    -> from the "R root" line (col2)

Input CSV must have a "Sample" column whose values match each report
file's sample name (e.g. sample "AI_1" -> "AI_1_bracken_report.txt").

Output: a new CSV with the same columns as the input, plus
classified_reads appended at the end.
"""

import csv
import os

# ---- CONFIG: edit these to match your setup ----
INPUT_CSV = "results/general_stats_1.csv"  # existing stats CSV to extend
BRACKEN_DIR = "bracken_results"  # folder containing *_bracken_report.txt files
FILE_PATTERN = "{sample}_bracken_report.txt"  # how sample name maps to report filename
OUTPUT_CSV = "results/general_stats_2.csv"
# ---------------------------------------------------


def parse_report(filepath):
    """Return classified_reads from a report file."""
    classified = 0

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue

            reads_in_clade = int(cols[1])
            rank_code = cols[3]

            if rank_code == "R":  # root = total classified reads
                classified = reads_in_clade
                break  # Stop once the root line is found

    return classified


def main():
    if not os.path.exists(INPUT_CSV):
        print(f"Input CSV not found: {INPUT_CSV}")
        return

    with open(INPUT_CSV, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    if "Sample" not in fieldnames:
        print(f"Input CSV has no 'Sample' column. Found columns: {fieldnames}")
        return

    out_fieldnames = fieldnames + ["classified_reads"]

    for row in rows:
        sample = row["Sample"].strip()
        report_path = os.path.join(
            BRACKEN_DIR, FILE_PATTERN.format(sample=sample)
        )

        if os.path.exists(report_path):
            classified = parse_report(report_path)
        else:
            print(
                f"WARNING: no report file found for sample '{sample}' (looked for {report_path}). "
                f"Leaving counts blank."
            )
            classified = ""

        row["classified_reads"] = classified

        if classified != "":
            print(f"{sample}: classified={classified:,}")

    os.makedirs(os.path.dirname(OUTPUT_CSV) or ".", exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nCSV saved to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
    
    
    