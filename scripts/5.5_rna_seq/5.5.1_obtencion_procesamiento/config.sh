#!/usr/bin/env bash
# =============================================================================
# config.sh — Central configuration for Fragaria x ananassa RNA-seq pipeline
# =============================================================================
set -euo pipefail

# ---- Project root (auto-detect from this script's location) -----------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Conda environment ------------------------------------------------------
ENV_NAME="rnaseq_aqp"

# ---- Hardware profile --------------------------------------------------------
# Detected: Intel i5-7400 (4 cores), 16 GB RAM
# Adjust these based on your system:
THREADS=4                   # Match physical cores (i5-7400 = 4 cores)
PARALLEL_JOBS=2             # Concurrent jobs for download / QC (limited by RAM)
SORT_MEMORY="1500M"         # Per-thread memory for samtools sort (conservative for 16GB)
SORT_THREADS=2              # Threads for samtools sort

# ---- Storage layout ----------------------------------------------------------
# SSD (440 GB, fast I/O): scripts, indices, references, results
# HDD (916 GB, large capacity): FASTQ, BAMs, trimmed reads, temp files
SSD_DIR="${PROJECT_DIR}"                          # /home/noe/work/RNA-seq_test (SSD)
# Override HDD_DIR via env var if your heavy data lives elsewhere:
#   export TFG_RNASEQ_HDD_DIR=/path/to/your/data
HDD_DIR="${TFG_RNASEQ_HDD_DIR:-/home/noe/work/datos/rnaseq_data}"

# ---- Reference files ---------------------------------------------------------
GENOME_GZ="${SSD_DIR}/Fragaria_ananassa_Benihoppe_Genome.fa.gz"
GFF3_GZ="${SSD_DIR}/Fragaria_ananassa_Benihoppe_Genome.gff3.gz"
GENOME="${SSD_DIR}/refs/genome.fa"
GFF3="${SSD_DIR}/refs/genome.gff3"
GTF="${SSD_DIR}/refs/genome.gtf"
TRANSCRIPTOME="${SSD_DIR}/refs/transcripts.fa"

# ---- Indices (SSD for fast access) -------------------------------------------
HISAT2_IDX="${SSD_DIR}/idx/genome"
SALMON_IDX="${SSD_DIR}/idx/salmon"

# ---- Heavy data directories (HDD) -------------------------------------------
FASTQ_DIR="${HDD_DIR}/fastq"
TRIMMED_DIR="${HDD_DIR}/trimmed"
BAM_DIR="${HDD_DIR}/bam"
SRA_CACHE_DIR="${HDD_DIR}/sra_cache"
TMP_DIR="${HDD_DIR}/tmp"

# ---- Light data directories (SSD) -------------------------------------------
QC_DIR="${SSD_DIR}/qc"
COUNTS_DIR="${SSD_DIR}/counts"
RESULTS_DIR="${SSD_DIR}/results"
LOG_DIR="${SSD_DIR}/logs"
DESIGN_DIR="${SSD_DIR}/design"
SALMON_CHECK_DIR="${SSD_DIR}/salmon_check"

# ---- Sample info -------------------------------------------------------------
SAMPLES_TSV="${DESIGN_DIR}/samples.tsv"
AQUAPORIN_TABLE="${SSD_DIR}/tabla_Aquaporinas_traduccion.tabular"

# ---- Parameters --------------------------------------------------------------
HISAT2_EXTRA_ARGS="--dta"   # Extra HISAT2 flags (--dta for StringTie compat)

# Strandedness: set to "auto" for auto-detection, or 0/1/2 if known
# 0 = unstranded, 1 = stranded (sense), 2 = reversely stranded
STRAND="auto"

# featureCounts parameters
FC_FEATURE_TYPE="exon"
FC_ATTR_TYPE="gene_id"

# ---- Data source selection ---------------------------------------------------
# Set to "sra" to download FASTQ from SRA and align locally
# Set to "galaxy" to download pre-aligned BAMs from Galaxy
DATA_SOURCE="sra"

# ---- Outlier flag ------------------------------------------------------------
# SRR30146487 (Roots Control rep 2) is a suspected outlier.
# Set to "include" to keep it in analysis, "exclude" to remove it.
OUTLIER_HANDLING="exclude"
OUTLIER_SRR="SRR30146487"
OUTLIER_SAMPLE="RootsCtrl_2"

# ---- Helper function: log with timestamp -------------------------------------
log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# ---- Helper function: check if conda env is active --------------------------
check_env() {
    if [[ "${CONDA_DEFAULT_ENV:-}" != "${ENV_NAME}" ]]; then
        echo "ERROR: Activate the environment first: micromamba activate ${ENV_NAME}"
        exit 1
    fi
}

# ---- Helper function: create all required directories ------------------------
ensure_dirs() {
    mkdir -p "${FASTQ_DIR}" "${TRIMMED_DIR}" "${BAM_DIR}" "${SRA_CACHE_DIR}" \
             "${TMP_DIR}" "${QC_DIR}" "${COUNTS_DIR}" "${RESULTS_DIR}" \
             "${LOG_DIR}" "${DESIGN_DIR}" "${SALMON_CHECK_DIR}" \
             "${SSD_DIR}/refs" "${SSD_DIR}/idx"
}
