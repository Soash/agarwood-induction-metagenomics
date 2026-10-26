#!/usr/bin/env python3
"""
fungi_AI_vs_NI_venn.py

Venn diagram of fungal species detected in:
    Artificially Induced (AI): AI_1, AI_2, AI_3  (union of species present in any)
    Naturally Induced   (NI): KNI, NI_1, NI_2   (union of species present in any)

A species counts as "present" in a group if it was detected (count_clade > 0)
in at least one sample belonging to that group. Control is excluded entirely,
same as fungi_AI_vs_NI_stats.py.

Uses the same Fungi-subtree parsing and EXCLUDED_SPECIES contaminant list as
fungi_summary.py / fungi_AI_vs_NI_stats.py, so results stay consistent across
all three scripts.

Outputs (written to OUTPUT_DIR):
    fungi_AI_vs_NI_venn.png          the Venn diagram
    fungi_AI_vs_NI_venn_species.csv  which species fall in each region (AI-only,
                                      NI-only, shared), one row per species

Run from the project root, e.g.:
    D:/geb/researches/Agarwood/final
    python scripts/fungi_AI_vs_NI_venn.py
"""

import csv
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
except ImportError:
    plt = None
    venn2 = None

# ============================================================
# CONFIG — adjust these if your layout differs
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # assumes script lives in scripts/

BRACKEN_REPORTS_DIR = PROJECT_ROOT / "bracken_results"
OUTPUT_DIR = PROJECT_ROOT / "fungi_summary"

VENN_OUT = OUTPUT_DIR / "fungi_AI_vs_NI_venn.png"
SPECIES_CSV_OUT = OUTPUT_DIR / "fungi_AI_vs_NI_venn_species.csv"

TAXID_FUNGI = 4751
TAXID_ROOT = 1

AI_GROUP = ["AI_1", "AI_2", "AI_3"]     # Artificially Induced
NI_GROUP = ["KNI", "NI_1", "NI_2"]      # Naturally Induced

EXCLUDED_SPECIES = {                   # species names to drop entirely from the output
    "Malassezia globosa",              # common human skin fungus, likely handler/lab contaminant
    "Malassezia restricta",            # common human skin fungus, likely handler/lab contaminant
    "Malassezia japonica",             # common human skin fungus, likely handler/lab contaminant
    "uncultured fungus",               # unresolved/ambiguous classification
    "uncultured Ascomycota",           # unresolved/ambiguous classification
    "fungal sp.",                      # unresolved/ambiguous classification
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
# STEP 2: Build the AI and NI species sets
# ============================================================

def build_group_species_set(group, all_species):
    species_set = set()
    for sample in group:
        for species, count in all_species.get(sample, {}).items():
            if count > 0:
                species_set.add(species)
    return species_set


# ============================================================
# STEP 3: Draw the Venn diagram
# ============================================================

def plot_venn(ai_species, ni_species):
    if plt is None or venn2 is None:
        print("  WARNING: matplotlib / matplotlib-venn not installed, skipping Venn diagram.")
        print("  Install with: pip install matplotlib matplotlib-venn")
        return

    fig, ax = plt.subplots(figsize=(7, 7))
    venn2(
        [ai_species, ni_species],
        set_labels=("Artificially Induced\n(AI_1, AI_2, AI_3)", "Naturally Induced\n(KNI, NI_1, NI_2)"),
        set_colors=("#4C72B0", "#DD8452"),
        alpha=0.6,
        ax=ax,
    )
    ax.set_title("Fungal Species Detected: AI vs NI")
    fig.savefig(VENN_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Venn diagram saved to: {VENN_OUT}")


# ============================================================
# STEP 4: Write CSV of which species fall in which region
# ============================================================

def write_species_csv(ai_species, ni_species):
    ai_only = sorted(ai_species - ni_species)
    ni_only = sorted(ni_species - ai_species)
    shared = sorted(ai_species & ni_species)

    with open(SPECIES_CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "region"])
        for sp in ai_only:
            writer.writerow([sp, "AI only"])
        for sp in ni_only:
            writer.writerow([sp, "NI only"])
        for sp in shared:
            writer.writerow([sp, "Shared (AI and NI)"])

    print(f"  Species-by-region CSV saved to: {SPECIES_CSV_OUT}")
    return ai_only, ni_only, shared


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
    for sample in all_samples:
        report_path = BRACKEN_REPORTS_DIR / f"{sample}_bracken_report.txt"
        if not report_path.exists():
            print(f"  WARNING: missing {report_path}")
            continue
        _, species_counts = parse_fungi_species(report_path)
        all_species[sample] = species_counts
        print(f"  {sample}: {len(species_counts)} fungal species")
    print()

    print("=" * 60)
    print("STEP 2: Building AI and NI species sets")
    print("=" * 60)
    ai_species = build_group_species_set(AI_GROUP, all_species)
    ni_species = build_group_species_set(NI_GROUP, all_species)
    print(f"  AI (union of AI_1, AI_2, AI_3): {len(ai_species)} distinct species")
    print(f"  NI (union of KNI, NI_1, NI_2):  {len(ni_species)} distinct species")
    print()

    print("=" * 60)
    print("STEP 3: Drawing Venn diagram")
    print("=" * 60)
    plot_venn(ai_species, ni_species)
    print()

    print("=" * 60)
    print("STEP 4: Writing species-by-region CSV")
    print("=" * 60)
    ai_only, ni_only, shared = write_species_csv(ai_species, ni_species)
    print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  AI only:      {len(ai_only)} species")
    print(f"  NI only:      {len(ni_only)} species")
    print(f"  Shared:       {len(shared)} species")
    print(f"  Total unique: {len(ai_species | ni_species)} species")


if __name__ == "__main__":
    main()
    
    
    