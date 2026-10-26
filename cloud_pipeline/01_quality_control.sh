#!/usr/bin/env bash
# ==============================================================================
# Step 01: Raw Read Quality Control and Preprocessing
# Infrastructure: Google Cloud Platform (GCP) Ubuntu Compute Instance
# ==============================================================================
set -euo pipefail

# 1. Directory Structure
mkdir -p raw_data fastqc_reports trimmed_data fastp_reports

# 2. Initial Quality Assessment with FastQC
echo "[INFO] Running FastQC on raw sequencing reads..."
fastqc raw_data/*.fastq.gz -o fastqc_reports -t 12

# 3. MultiQC Aggregation of Raw Reads
echo "[INFO] Generating MultiQC summary for raw reads..."
multiqc fastqc_reports -o fastqc_reports/multiqc_raw

# 4. Adapter Trimming and Poly-G Quality Filtering via fastp
echo "[INFO] Executing fastp quality filtering..."
for R1 in raw_data/*_R1.fastq.gz; do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    
    fastp \
        --in1 "$R1" --in2 "$R2" \
        --out1 "trimmed_data/${SAMPLE}_R1.trimmed.fastq.gz" \
        --out2 "trimmed_data/${SAMPLE}_R2.trimmed.fastq.gz" \
        --html "fastp_reports/${SAMPLE}_fastp.html" \
        --json "fastp_reports/${SAMPLE}_fastp.json" \
        --thread 8 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50
        
    echo "[DONE] Trimmed: $SAMPLE"
done

# 5. MultiQC on Trimmed Reads
echo "[INFO] Generating post-trimming MultiQC report..."
multiqc fastp_reports -o fastp_reports/multiqc_trimmed
echo "[SUCCESS] Quality control and trimming complete."
