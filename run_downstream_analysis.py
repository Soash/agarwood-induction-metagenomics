#!/usr/bin/env python3
"""
Master Downstream Analysis Orchestrator for Agarwood Metagenomics
Executes the downstream statistical, CAZyme, taxonomic, and visualization pipelines.
"""

import subprocess
import sys
import os

SCRIPTS_ORDER = [
    ("scripts/add_classified_counts.py", "Computing classified vs unclassified read metrics"),
    ("scripts/cazyme_pipeline.py", "Processing CAZyme profiles and CPM abundance tables"),
    ("scripts/domain_summary.py", "Generating multi-domain taxonomic composition metrics"),
    ("scripts/fungi_summary.py", "Computing fungal species counts, richness, and CPM abundances"),
    ("scripts/fungi_AI_vs_NI_stats.py", "Calculating differential abundance statistics (AI vs NI)"),
    ("scripts/plot_agarwood_predictions.py", "Generating AWI-Pred decision tree and CAZyme scatter plots"),
    ("scripts/fungi_AI_vs_NI_venn.py", "Constructing species overlap Venn diagrams (AI vs NI)"),
    ("scripts/fungi_species_heatmap.py", "Rendering fungal species abundance heatmaps")
]

def run_script(script_path, description):
    print(f"\n{'='*70}")
    print(f"[RUNNING] {script_path}: {description}")
    print(f"{'='*70}")
    
    if not os.path.exists(script_path):
        print(f"[ERROR] Script not found: {script_path}", file=sys.stderr)
        return False
        
    res = subprocess.run([sys.executable, script_path], text=True)
    if res.returncode != 0:
        print(f"[FAILED] {script_path} exited with code {res.returncode}", file=sys.stderr)
        return False
    print(f"[SUCCESS] {script_path} completed successfully.")
    return True

def main():
    print("Starting Agarwood Metagenomics Downstream Analysis Pipeline...")
    for script_path, desc in SCRIPTS_ORDER:
        success = run_script(script_path, desc)
        if not success:
            print(f"\n[ABORT] Pipeline halted due to error in {script_path}")
            sys.exit(1)
            
    print("\n" + "="*70)
    print("[COMPLETED] All downstream analyses and figures generated successfully.")
    print("="*70)

if __name__ == "__main__":
    main()
