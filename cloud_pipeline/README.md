# Google Cloud Platform (GCP) Ubuntu Compute Execution Guide

This directory contains production shell scripts utilized for executing heavy metagenomic and computational workloads on a Google Cloud Platform Compute Engine Ubuntu instance.

## Cloud Environment Specifications
* **Operating System:** Ubuntu 22.04 LTS (x86_64)
* **vCPUs / RAM:** 16 vCPUs / 64 GB RAM
* **Persistent Disk:** 500 GB SSD mounted for raw fastq processing; external HDD mounted at `/mnt/disks/hdd/` hosting the uncompressed Kraken2/Bracken reference database (`k2_core_nt`, ~100 GB).

## Sequential Execution Pipeline
1. **`01_quality_control.sh`**: Executes FastQC on raw reads, runs MultiQC, and applies `fastp` adapter clipping and poly-G filtering ($Q \ge 20$).
2. **`02_host_depletion.sh`**: Downloads canonical *Aquilaria malaccensis* and human contaminant reference genomes, builds a combined Bowtie2 index, maps trimmed reads, and isolates unmapped non-host microbial pairs.
3. **`03_taxonomic_profiling.sh`**: Classifies microbial reads across domains using `Kraken2` and re-estimates species-level relative abundance using `Bracken v3.0.1`.
4. **`04_dbcan_cazyme.sh`**: Executes HMMER and DIAMOND searches against the CAZy database via `dbCAN` to quantify Glycoside Hydrolases (GH), Auxiliary Activities (AA), and carbohydrate-binding modules.
