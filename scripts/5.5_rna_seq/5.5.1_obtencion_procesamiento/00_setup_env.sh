#!/usr/bin/env bash
# =============================================================================
# 00_setup_env.sh — Create micromamba environment for RNA-seq pipeline
# =============================================================================
set -euo pipefail

ENV_NAME="rnaseq_aqp"

echo "=== Creating micromamba environment: ${ENV_NAME} ==="

# Remove existing env if present
micromamba env list | grep -q "${ENV_NAME}" && micromamba env remove -n "${ENV_NAME}" -y || true

# Create environment with all required tools
micromamba create -n "${ENV_NAME}" -c conda-forge -c bioconda -y \
    sra-tools=3.1 \
    fastp \
    hisat2 \
    samtools \
    subread \
    salmon \
    gffread \
    parallel \
    pigz \
    multiqc \
    r-base \
    r-tidyverse \
    r-pheatmap \
    r-rcolorbrewer \
    r-ggrepel \
    r-gplots \
    r-optparse \
    bioconductor-deseq2 \
    bioconductor-edger \
    bioconductor-tximport \
    bioconductor-enhancedvolcano \
    python=3.12 \
    pandas \
    matplotlib \
    seaborn

echo "=== Environment '${ENV_NAME}' created successfully ==="
echo "Activate with: micromamba activate ${ENV_NAME}"
