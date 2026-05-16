#!/usr/bin/env bash
# =============================================================================
# 01_prepare_refs.sh — Decompress references, convert GFF3→GTF, build indices
#
# HARDWARE NOTES (i5-7400, 16 GB RAM):
#   - HISAT2-build with splice sites on ~793 MB genome needs ~20+ GB RAM
#   - On 16 GB systems, we create temporary swap on HDD first
#   - In galaxy mode, HISAT2 index is SKIPPED (BAMs are pre-aligned)
#   - Salmon index is always built (needed for strandedness check)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

# Read DATA_MODE from environment (set by run_pipeline.sh)
DATA_MODE="${DATA_MODE:-sra}"

log_msg "=== Step 1: Preparing reference files ==="
log_msg "Data mode: ${DATA_MODE}"

# ---- 1.1 Decompress genome --------------------------------------------------
if [[ ! -f "${GENOME}" ]]; then
    log_msg "Decompressing genome..."
    gunzip -ck "${GENOME_GZ}" > "${GENOME}"
    samtools faidx "${GENOME}"
    log_msg "Genome decompressed: $(wc -c < "${GENOME}") bytes"
else
    log_msg "Genome already exists, skipping decompression."
fi

# ---- 1.2 Decompress GFF3 ----------------------------------------------------
if [[ ! -f "${GFF3}" ]]; then
    log_msg "Decompressing GFF3..."
    gunzip -ck "${GFF3_GZ}" > "${GFF3}"
    log_msg "GFF3 decompressed."
else
    log_msg "GFF3 already exists, skipping."
fi

# ---- 1.3 Convert GFF3 to GTF ------------------------------------------------
if [[ ! -f "${GTF}" ]]; then
    log_msg "Converting GFF3 to GTF with gffread..."
    gffread "${GFF3}" -T -o "${GTF}"

    # Verify GTF has gene_id attribute
    GENE_ID_COUNT=$(grep -c 'gene_id' "${GTF}" || true)
    if [[ "${GENE_ID_COUNT}" -eq 0 ]]; then
        log_msg "ERROR: GTF conversion failed — no gene_id attributes found."
        exit 1
    fi
    log_msg "GTF created with ${GENE_ID_COUNT} lines containing gene_id."
else
    log_msg "GTF already exists, skipping."
fi

# ---- 1.4 Extract transcriptome for salmon -----------------------------------
if [[ ! -f "${TRANSCRIPTOME}" ]]; then
    log_msg "Extracting transcriptome sequences with gffread..."
    gffread "${GFF3}" -g "${GENOME}" -w "${TRANSCRIPTOME}"
    log_msg "Transcriptome: $(grep -c '^>' "${TRANSCRIPTOME}") sequences."
else
    log_msg "Transcriptome already exists, skipping."
fi

# ---- 1.5 Build HISAT2 index -------------------------------------------------
# SKIP in galaxy mode — pre-aligned BAMs don't need local HISAT2 index
if [[ "${DATA_MODE}" == "galaxy" ]]; then
    log_msg "Galaxy mode: SKIPPING HISAT2 index (BAMs are pre-aligned)."
elif [[ -f "${HISAT2_IDX}.1.ht2" ]] || [[ -f "${HISAT2_IDX}.1.ht2l" ]]; then
    log_msg "HISAT2 index already exists, skipping."
else
    log_msg "Extracting splice sites and exons from GTF..."
    hisat2_extract_splice_sites.py "${GTF}" > "${SSD_DIR}/refs/splicesites.tsv"
    hisat2_extract_exons.py "${GTF}" > "${SSD_DIR}/refs/exons.tsv"

    # Check available RAM
    AVAIL_MEM_GB=$(awk '/MemAvailable/ {printf "%.0f", $2/1024/1024}' /proc/meminfo)
    log_msg "Available RAM: ${AVAIL_MEM_GB} GB"

    # The ~793 MB octoploid genome with splice sites needs ~20 GB+ RAM for HISAT2-build.
    # On 16 GB systems, we MUST add temporary swap space.
    SWAP_FILE="${HDD_DIR}/tmp/hisat2_swap"
    SWAP_CREATED=false

    if [[ "${AVAIL_MEM_GB}" -lt 20 ]]; then
        log_msg "WARNING: Only ${AVAIL_MEM_GB} GB available, need ~20 GB for HISAT2-build."
        log_msg "Creating temporary 20 GB swap file on HDD..."

        if [[ ! -f "${SWAP_FILE}" ]]; then
            # Create swap file on HDD (large, slow, but better than OOM)
            sudo dd if=/dev/zero of="${SWAP_FILE}" bs=1M count=20480 status=progress 2>&1 || {
                log_msg "ERROR: Cannot create swap file. Run with sudo or create swap manually:"
                log_msg "  sudo dd if=/dev/zero of=${SWAP_FILE} bs=1M count=20480"
                log_msg "  sudo chmod 600 ${SWAP_FILE}"
                log_msg "  sudo mkswap ${SWAP_FILE}"
                log_msg "  sudo swapon ${SWAP_FILE}"
                log_msg "Then re-run the pipeline."
                exit 1
            }
            sudo chmod 600 "${SWAP_FILE}"
            sudo mkswap "${SWAP_FILE}" 2>&1
            sudo swapon "${SWAP_FILE}" 2>&1
            SWAP_CREATED=true
            log_msg "Swap file activated. Total virtual memory now:"
            free -h | head -3
        elif swapon --show | grep -q "${SWAP_FILE}"; then
            log_msg "Swap file already active."
        else
            sudo swapon "${SWAP_FILE}" 2>&1 || {
                log_msg "WARNING: Could not activate existing swap file."
                log_msg "Creating a new swap file..."
                sudo dd if=/dev/zero of="${SWAP_FILE}" bs=1M count=20480 status=progress 2>&1
                sudo chmod 600 "${SWAP_FILE}"
                sudo mkswap "${SWAP_FILE}" 2>&1
                sudo swapon "${SWAP_FILE}" 2>&1
                SWAP_CREATED=true
            }
        fi
    fi

    # Build with 1 thread to minimize RAM overhead
    INDEX_THREADS=1
    log_msg "Building HISAT2 index with ${INDEX_THREADS} thread (RAM-safe + swap)..."
    log_msg "This may take 60-120 minutes with swap on HDD. Be patient..."
    hisat2-build \
        --ss "${SSD_DIR}/refs/splicesites.tsv" \
        --exon "${SSD_DIR}/refs/exons.tsv" \
        -p "${INDEX_THREADS}" \
        "${GENOME}" \
        "${HISAT2_IDX}" \
        2>&1 | tee "${LOG_DIR}/hisat2_build.log"
    log_msg "HISAT2 index built."

    # Clean up temporary swap
    if [[ "${SWAP_CREATED}" == true ]]; then
        log_msg "Removing temporary swap file..."
        sudo swapoff "${SWAP_FILE}" 2>/dev/null || true
        rm -f "${SWAP_FILE}" 2>/dev/null || true
        log_msg "Swap file removed."
    fi
fi

# ---- 1.6 Build salmon index -------------------------------------------------
if [[ ! -d "${SALMON_IDX}" ]]; then
    log_msg "Building salmon index..."
    salmon index \
        -t "${TRANSCRIPTOME}" \
        -i "${SALMON_IDX}" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/salmon_index.log"
    log_msg "Salmon index built."
else
    log_msg "Salmon index already exists, skipping."
fi

log_msg "=== Step 1 complete ==="
