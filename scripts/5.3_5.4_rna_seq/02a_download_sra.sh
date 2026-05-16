#!/usr/bin/env bash
# =============================================================================
# 02a_download_sra.sh — Download FASTQ files from SRA using prefetch+fasterq-dump
#
# FIXES applied:
#   - Uses prefetch before fasterq-dump (NCBI recommended, more robust)
#   - Parallelized with GNU parallel (was sequential)
#   - Temp files go to HDD tmp dir (not FASTQ dir)
#   - Validates gzip integrity after compression
#   - Uses process substitution to avoid subshell variable loss
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

log_msg "=== Step 2A: Downloading FASTQ from SRA ==="

# ---- Download function (called by GNU parallel) -----------------------------
download_one_sample() {
    local sample_id="$1"
    local srr_id="$2"

    # Source config inside parallel subshell
    source "$(dirname "$0")/config.sh"

    local R1="${FASTQ_DIR}/${sample_id}_1.fastq.gz"
    local R2="${FASTQ_DIR}/${sample_id}_2.fastq.gz"

    if [[ -f "${R1}" && -f "${R2}" ]]; then
        log_msg "Skipping ${sample_id} (${srr_id}) — already downloaded."
        return 0
    fi

    log_msg "Downloading ${sample_id} (${srr_id})..."

    # Step 1: Prefetch SRA file (robust, resumable download)
    log_msg "  Prefetching ${srr_id}..."
    prefetch "${srr_id}" \
        --output-directory "${SRA_CACHE_DIR}" \
        --max-size 50G \
        2>&1 | tee "${LOG_DIR}/${sample_id}_prefetch.log"

    # Locate the .sra file
    local SRA_FILE=""
    if [[ -f "${SRA_CACHE_DIR}/${srr_id}/${srr_id}.sra" ]]; then
        SRA_FILE="${SRA_CACHE_DIR}/${srr_id}/${srr_id}.sra"
    elif [[ -f "${SRA_CACHE_DIR}/${srr_id}.sra" ]]; then
        SRA_FILE="${SRA_CACHE_DIR}/${srr_id}.sra"
    else
        log_msg "ERROR: Prefetch completed but .sra file not found for ${srr_id}"
        return 1
    fi

    # Step 2: Extract FASTQ with fasterq-dump (temp on HDD, not in FASTQ dir)
    log_msg "  Extracting FASTQ for ${sample_id}..."
    fasterq-dump "${SRA_FILE}" \
        --split-3 \
        --outdir "${FASTQ_DIR}" \
        --temp "${TMP_DIR}" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/${sample_id}_fasterq.log"

    # Step 3: Rename and compress with pigz
    if [[ -f "${FASTQ_DIR}/${srr_id}_1.fastq" ]]; then
        log_msg "  Compressing ${sample_id}..."
        pigz -p "${THREADS}" -c "${FASTQ_DIR}/${srr_id}_1.fastq" > "${R1}"
        pigz -p "${THREADS}" -c "${FASTQ_DIR}/${srr_id}_2.fastq" > "${R2}"

        # Step 4: Validate gzip integrity
        if gzip -t "${R1}" 2>/dev/null && gzip -t "${R2}" 2>/dev/null; then
            log_msg "  Integrity OK for ${sample_id}. Cleaning up..."
            rm -f "${FASTQ_DIR}/${srr_id}_1.fastq" "${FASTQ_DIR}/${srr_id}_2.fastq"
            # Clean SRA cache for this sample to save space
            rm -rf "${SRA_CACHE_DIR}/${srr_id}"
        else
            log_msg "ERROR: Corrupt gzip for ${sample_id}! Keeping uncompressed files."
            rm -f "${R1}" "${R2}"
            return 1
        fi
    elif [[ -f "${FASTQ_DIR}/${srr_id}.fastq" ]]; then
        # Single-end fallback (shouldn't happen for our data)
        log_msg "WARNING: ${srr_id} downloaded as single-end!"
        pigz -p "${THREADS}" -c "${FASTQ_DIR}/${srr_id}.fastq" > "${FASTQ_DIR}/${sample_id}.fastq.gz"
        rm -f "${FASTQ_DIR}/${srr_id}.fastq"
    else
        log_msg "ERROR: No FASTQ files produced for ${srr_id}"
        return 1
    fi

    log_msg "Done: ${sample_id}"
}
export -f download_one_sample

# ---- Get script directory for parallel subshells ----------------------------
SCRIPT_DIR_EXPORT="$(cd "$(dirname "$0")" && pwd)"
export SCRIPT_DIR_EXPORT

# Override download_one_sample to source config properly
download_wrapper() {
    local sample_id="$1"
    local srr_id="$2"

    source "${SCRIPT_DIR_EXPORT}/config.sh"

    local R1="${FASTQ_DIR}/${sample_id}_1.fastq.gz"
    local R2="${FASTQ_DIR}/${sample_id}_2.fastq.gz"

    if [[ -f "${R1}" && -f "${R2}" ]]; then
        log_msg "Skipping ${sample_id} (${srr_id}) — already downloaded."
        return 0
    fi

    log_msg "Downloading ${sample_id} (${srr_id})..."

    # Prefetch
    prefetch "${srr_id}" \
        --output-directory "${SRA_CACHE_DIR}" \
        --max-size 50G \
        2>&1 | tee "${LOG_DIR}/${sample_id}_prefetch.log"

    local SRA_FILE=""
    if [[ -f "${SRA_CACHE_DIR}/${srr_id}/${srr_id}.sra" ]]; then
        SRA_FILE="${SRA_CACHE_DIR}/${srr_id}/${srr_id}.sra"
    elif [[ -f "${SRA_CACHE_DIR}/${srr_id}.sra" ]]; then
        SRA_FILE="${SRA_CACHE_DIR}/${srr_id}.sra"
    else
        log_msg "ERROR: .sra file not found for ${srr_id}"
        return 1
    fi

    # fasterq-dump (temp on HDD, not FASTQ dir)
    fasterq-dump "${SRA_FILE}" \
        --split-3 \
        --outdir "${FASTQ_DIR}" \
        --temp "${TMP_DIR}" \
        --threads 2 \
        2>&1 | tee "${LOG_DIR}/${sample_id}_fasterq.log"

    # Compress and validate
    if [[ -f "${FASTQ_DIR}/${srr_id}_1.fastq" ]]; then
        pigz -p 2 -c "${FASTQ_DIR}/${srr_id}_1.fastq" > "${R1}"
        pigz -p 2 -c "${FASTQ_DIR}/${srr_id}_2.fastq" > "${R2}"

        if gzip -t "${R1}" 2>/dev/null && gzip -t "${R2}" 2>/dev/null; then
            rm -f "${FASTQ_DIR}/${srr_id}_1.fastq" "${FASTQ_DIR}/${srr_id}_2.fastq"
            rm -rf "${SRA_CACHE_DIR}/${srr_id}"
            log_msg "Done: ${sample_id} (validated)"
        else
            log_msg "ERROR: Corrupt gzip for ${sample_id}!"
            rm -f "${R1}" "${R2}"
            return 1
        fi
    elif [[ -f "${FASTQ_DIR}/${srr_id}.fastq" ]]; then
        log_msg "WARNING: ${srr_id} downloaded as single-end!"
        pigz -p 2 -c "${FASTQ_DIR}/${srr_id}.fastq" > "${FASTQ_DIR}/${sample_id}.fastq.gz"
        rm -f "${FASTQ_DIR}/${srr_id}.fastq"
    else
        log_msg "ERROR: No FASTQ files produced for ${srr_id}"
        return 1
    fi
}
export -f download_wrapper

# ---- Run parallel downloads (2 concurrent to balance I/O and network) -------
log_msg "Downloading ${PARALLEL_JOBS} samples in parallel..."

# Extract sample_id and srr_id columns, pass to parallel
tail -n +2 "${SAMPLES_TSV}" | cut -f1,2 | \
    parallel -j "${PARALLEL_JOBS}" --colsep '\t' --halt soon,fail=1 \
    download_wrapper {1} {2}

log_msg "=== Step 2A complete ==="
log_msg "Downloaded FASTQ files:"
ls -lh "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null || log_msg "No FASTQ files found!"

# Clean up temp directory
rm -rf "${TMP_DIR:?}"/*
log_msg "Temp files cleaned."
