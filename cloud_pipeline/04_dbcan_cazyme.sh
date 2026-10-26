#!/usr/bin/env bash
# ==============================================================================
# Step 04: dbCAN Carbohydrate-Active Enzymes (CAZyme) Profiling
# Infrastructure: Google Cloud Platform (GCP) Ubuntu Compute Instance
# ==============================================================================
set -euo pipefail

mkdir -p dbCAN/db dbCAN/results/class_counts dbCAN/results/filtered

# 1. Download and Format dbCAN Database
echo "[INFO] Initializing dbCAN database..."
# dbcan_build --cpus 8 --db-dir dbCAN/db --clean

# 2. Screening Fungal Proteomes / Metagenomic Contigs Against CAZy
echo "[INFO] Running run_dbcan against CAZy database..."
for PROTEIN_FAA in data/fungal_proteomes/*.faa; do
    [ -e "$PROTEIN_FAA" ] || continue
    SP_NAME=$(basename "$PROTEIN_FAA" .faa)
    
    # Sanitize fasta headers
    sed 's/:/_/g' "$PROTEIN_FAA" > "${PROTEIN_FAA}.clean"
    
    run_dbcan "${PROTEIN_FAA}.clean" protein \
        --out_dir "dbCAN/results/${SP_NAME}_dbcan" \
        --db_dir dbCAN/db \
        --cpu 8
        
    rm -f "${PROTEIN_FAA}.clean"
    echo "[DONE] dbCAN profiled: $SP_NAME"
done

echo "[SUCCESS] dbCAN CAZyme profiling complete."
