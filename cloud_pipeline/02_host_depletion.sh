#!/usr/bin/env bash
# ==============================================================================
# Step 02: Host Genome Depletion (Aquilaria malaccensis + Human Contaminant)
# Infrastructure: Google Cloud Platform (GCP) Ubuntu Compute Instance
# ==============================================================================
set -euo pipefail

mkdir -p host_genome host_removed

# 1. Download Reference Genomes using NCBI Datasets CLI
echo "[INFO] Downloading Aquilaria malaccensis reference genome (GCA_037177235.1)..."
datasets download genome accession GCA_037177235.1 --include genome --filename host_genome/aquilaria.zip
unzip -q host_genome/aquilaria.zip -d host_genome/aquilaria
mv host_genome/aquilaria/ncbi_dataset/data/GCA_037177235.1/*.fna host_genome/Aquilaria_malaccensis.fna
rm -rf host_genome/aquilaria host_genome/aquilaria.zip

echo "[INFO] Downloading Homo sapiens contaminant reference (GCF_000001405.40)..."
datasets download genome accession GCF_000001405.40 --include genome --filename host_genome/human.zip
unzip -q host_genome/human.zip -d host_genome/human
mv host_genome/human/ncbi_dataset/data/GCF_000001405.40/*.fna host_genome/Homo_sapiens.fna
rm -rf host_genome/human host_genome/human.zip

# 2. Build Combined Dual-Host Bowtie2 Index
echo "[INFO] Concatenating reference genomes and indexing with Bowtie2..."
cat host_genome/Aquilaria_malaccensis.fna host_genome/Homo_sapiens.fna > host_genome/combined_host.fna
bowtie2-build --threads 12 host_genome/combined_host.fna host_genome/combined_host_index

# 3. Map Trimmed Reads and Extract Pure Microbial Unmapped Read Pairs
echo "[INFO] Aligning reads to combined host index and extracting unmapped microbial pairs..."
for R1 in trimmed_data/*_R1.trimmed.fastq.gz; do
    R2="${R1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}"
    SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)
    
    bowtie2 -x host_genome/combined_host_index \
        -1 "$R1" -2 "$R2" \
        --un-conc-gz "host_removed/${SAMPLE}_unmapped_%.fastq.gz" \
        -S /dev/null \
        -p 12 \
        --very-sensitive \
        2> "host_removed/${SAMPLE}_bowtie2.log"
        
    echo "[DONE] Host-depleted: $SAMPLE"
done

echo "[SUCCESS] Host depletion complete. Non-host reads ready for taxonomic classification."
