#!/usr/bin/env bash
# =============================================================================
# 04_align.sh — HISAT2 alignment of trimmed reads to genome
#
# FIXES applied:
#   - samtools sort uses TMP_DIR on HDD for temp files (not RAM-only)
#   - Reduced sort memory to fit 16GB system
#   - Explicit -T flag for samtools sort temp prefix
#   - Uses process substitution instead of pipe|while
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

log_msg "=== Step 4: HISAT2 alignment ==="

# Verify HISAT2 index exists
if [[ ! -f "${HISAT2_IDX}.1.ht2" ]] && [[ ! -f "${HISAT2_IDX}.1.ht2l" ]]; then
    log_msg "ERROR: HISAT2 index not found. Run 01_prepare_refs.sh first."
    exit 1
fi

# Read sample info using process substitution (avoids subshell variable loss)
while IFS=$'\t' read -r sample_id srr_id tissue condition replicate outlier_flag; do

    BAM_FILE="${BAM_DIR}/${sample_id}.bam"
    R1="${TRIMMED_DIR}/${sample_id}_1.fastq.gz"
    R2="${TRIMMED_DIR}/${sample_id}_2.fastq.gz"

    if [[ -f "${BAM_FILE}" ]]; then
        log_msg "Skipping ${sample_id} — BAM already exists."
        continue
    fi

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        log_msg "WARNING: Trimmed reads not found for ${sample_id}, skipping."
        continue
    fi

    log_msg "Aligning ${sample_id}..."

    # HISAT2 alignment piped to samtools sort
    # Memory: HISAT2 ~6-8GB for this genome + samtools sort ~3GB
    # Total ~11GB, fits in 16GB with headroom
    hisat2 \
        -x "${HISAT2_IDX}" \
        ${HISAT2_EXTRA_ARGS} \
        -1 "${R1}" \
        -2 "${R2}" \
        --threads "${THREADS}" \
        --summary-file "${LOG_DIR}/${sample_id}_hisat2.log" \
        --new-summary \
    | samtools sort \
        -@ "${SORT_THREADS}" \
        -m "${SORT_MEMORY}" \
        -T "${TMP_DIR}/sort_${sample_id}" \
        -o "${BAM_FILE}"

    samtools index "${BAM_FILE}"
    log_msg "Done: ${sample_id} — $(samtools flagstat "${BAM_FILE}" | head -1)"

done < <(tail -n +2 "${SAMPLES_TSV}")

# ---- Alignment rate summary --------------------------------------------------
log_msg "=== Alignment rate summary ==="
echo -e "sample\ttotal_reads\taligned\talignment_rate" > "${LOG_DIR}/alignment_summary.tsv"

for log_file in "${LOG_DIR}"/*_hisat2.log; do
    [[ -f "${log_file}" ]] || continue
    sample=$(basename "${log_file}" _hisat2.log)
    # Parse HISAT2 new-summary format
    total=$(grep "Total pairs" "${log_file}" | awk '{print $1}' || echo "N/A")
    rate=$(grep "Overall alignment rate" "${log_file}" | awk '{print $1}' || echo "N/A")
    aligned=$(grep "Aligned concordantly" "${log_file}" | head -1 | awk '{print $1}' || echo "N/A")
    echo -e "${sample}\t${total}\t${aligned}\t${rate}" >> "${LOG_DIR}/alignment_summary.tsv"
done

cat "${LOG_DIR}/alignment_summary.tsv" | column -t
log_msg "=== Step 4 complete ==="
