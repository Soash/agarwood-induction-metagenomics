#!/usr/bin/env bash
# ==============================================================================
# Step 03: Taxonomic Profiling with Kraken2 and Bracken
# Infrastructure: Google Cloud Platform (GCP) Ubuntu Compute Instance
# Reference DB: Kraken2 PlusPF / Core-NT Database mounted on persistent disk
# ==============================================================================
set -euo pipefail

DB_DIR="/mnt/disks/hdd/k2_core_nt"
mkdir -p kraken2_results bracken_results

# 1. Kraken2 Taxonomic Classification
echo "[INFO] Running Kraken2 classification against k2_core_nt database..."
for R1 in host_removed/*_unmapped_1.fastq.gz; do
    R2="${R1/_unmapped_1.fastq.gz/_unmapped_2.fastq.gz}"
    SAMPLE=$(basename "$R1" _unmapped_1.fastq.gz)
    
    echo "[INFO] Processing sample: $SAMPLE"
    kraken2 --db "$DB_DIR" \
        --threads 12 \
        --paired "$R1" "$R2" \
        --report "kraken2_results/${SAMPLE}_report.txt" \
        --output "kraken2_results/${SAMPLE}_output.txt"
done

# 2. Bracken Abundance Re-estimation at Species Level
echo "[INFO] Re-estimating species-level abundances with Bracken..."
declare -A READLEN
READLEN[AI_1]=150
READLEN[AI_2]=100
READLEN[AI_3]=150
READLEN[Control]=100
READLEN[KNI]=100
READLEN[NI_1]=100
READLEN[NI_2]=100

for report in kraken2_results/*_report.txt; do
    SAMPLE=$(basename "$report" _report.txt)
    RLEN=${READLEN[$SAMPLE]:-100}
    
    echo "[INFO] Running Bracken for $SAMPLE (read length: $RLEN)..."
    bracken -d "$DB_DIR" \
        -i "$report" \
        -o "bracken_results/${SAMPLE}_bracken.txt" \
        -w "bracken_results/${SAMPLE}_bracken_report.txt" \
        -r "$RLEN" \
        -l S \
        -t 10
done

echo "[SUCCESS] Taxonomic classification and abundance estimation complete."
