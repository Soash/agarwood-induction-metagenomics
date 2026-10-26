#!/usr/bin/env python3
"""
other_domains_species_richness.py

Calculates and plots species richness (number of distinct species detected)
for Bacteria, Archaea, and Viruses, similar to fungi_summary.py.
"""

import csv
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
except ImportError:
    plt = None

PROJECT_ROOT = Path(__file__).resolve().parent.parent

BRACKEN_REPORTS_DIR = PROJECT_ROOT / "bracken_results"
OUTPUT_DIR = PROJECT_ROOT / "domain_summary"

TAXID_BACTERIA = 2
TAXID_ARCHAEA = 2157
TAXID_VIRUSES = 10239

DOMAINS = {
    "Bacteria": TAXID_BACTERIA,
    "Archaea": TAXID_ARCHAEA,
    "Virus": TAXID_VIRUSES
}

def parse_domain_species(path: Path, target_taxid: int):
    """Return {species_name: count_clade} for a specific domain."""
    species_counts = {}

    within_domain = False
    domain_indent = None

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

            if not within_domain:
                if tax_id == target_taxid:
                    within_domain = True
                    domain_indent = indent
                continue

            # Already inside the target domain subtree
            if indent <= domain_indent:
                within_domain = False  # exited subtree; keep scanning
                continue

            if rank_code == "S":
                species_counts[name] = species_counts.get(name, 0) + count_clade

    return species_counts

def plot_species_richness_by_sample(all_species, samples, domain_name):
    if plt is None:
        print(f"  WARNING: matplotlib not installed, skipping {domain_name} species richness chart.")
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

    ax.set_ylabel(f"Number of distinct {domain_name} species detected")
    ax.set_xlabel("Sample")
    ax.set_title(f"{domain_name} Species Richness by Sample")
    
    out_path = OUTPUT_DIR / f"{domain_name.lower()}_species_richness_by_sample.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Species richness bar chart saved to: {out_path}")

def write_richness_csv(all_domain_species, samples):
    out_path = OUTPUT_DIR / "other_domains_species_richness.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Sample"] + list(DOMAINS.keys())
        writer.writerow(header)
        
        for sample in samples:
            row = [sample]
            for domain_name in DOMAINS.keys():
                richness = len(all_domain_species.get(domain_name, {}).get(sample, {}))
                row.append(richness)
            writer.writerow(row)
    print(f"Species richness CSV saved to: {out_path}")

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    samples = sorted(
        p.name.replace("_bracken_report.txt", "")
        for p in BRACKEN_REPORTS_DIR.glob("*_bracken_report.txt")
    )
    print(f"Detected samples: {samples}\n")

    all_domain_species = {}

    for domain_name, tax_id in DOMAINS.items():
        print("=" * 60)
        print(f"Parsing species for {domain_name} (TaxID: {tax_id})")
        print("=" * 60)
        
        all_species = {}
        for sample in samples:
            report_path = BRACKEN_REPORTS_DIR / f"{sample}_bracken_report.txt"
            if not report_path.exists():
                print(f"  WARNING: missing {report_path}")
                continue
            
            species_counts = parse_domain_species(report_path, tax_id)
            all_species[sample] = species_counts
            print(f"  {sample}: {len(species_counts)} {domain_name} species")
        
        all_domain_species[domain_name] = all_species
        plot_species_richness_by_sample(all_species, samples, domain_name)
        print()
        
    print("=" * 60)
    print("Writing CSV summary")
    print("=" * 60)
    write_richness_csv(all_domain_species, samples)

if __name__ == "__main__":
    main()
