#!/usr/bin/env bash
# =============================================================================
# 03_qc_trim.sh — Quality control and adapter trimming with fastp
#
# FIXES applied:
#   - Parallelized with GNU parallel (2 samples at a time for 4-core CPU)
#   - Uses process substitution instead of pipe|while (avoid subshell)
#   - fastp threads reduced to 2 per sample when running 2 in parallel
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

log_msg "=== Step 3: QC and trimming with fastp ==="

SCRIPT_DIR_QC="$(cd "$(dirname "$0")" && pwd)"
export SCRIPT_DIR_QC

# ---- Trim function (called by GNU parallel) ---------------------------------
trim_one_sample() {
    local sample_id="$1"
    source "${SCRIPT_DIR_QC}/config.sh"

    local R1_IN="${FASTQ_DIR}/${sample_id}_1.fastq.gz"
    local R2_IN="${FASTQ_DIR}/${sample_id}_2.fastq.gz"
    local R1_OUT="${TRIMMED_DIR}/${sample_id}_1.fastq.gz"
    local R2_OUT="${TRIMMED_DIR}/${sample_id}_2.fastq.gz"
    local HTML_REPORT="${QC_DIR}/${sample_id}_fastp.html"
    local JSON_REPORT="${QC_DIR}/${sample_id}_fastp.json"

    if [[ -f "${R1_OUT}" && -f "${R2_OUT}" ]]; then
        log_msg "Skipping ${sample_id} — already trimmed."
        return 0
    fi

    if [[ ! -f "${R1_IN}" || ! -f "${R2_IN}" ]]; then
        log_msg "WARNING: FASTQ files not found for ${sample_id}, skipping."
        return 0
    fi

    log_msg "Running fastp on ${sample_id}..."
    fastp \
        -i "${R1_IN}" -I "${R2_IN}" \
        -o "${R1_OUT}" -O "${R2_OUT}" \
        --html "${HTML_REPORT}" \
        --json "${JSON_REPORT}" \
        --thread 2 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36 \
        --cut_front --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        2>&1 | tee "${LOG_DIR}/${sample_id}_fastp.log"

    log_msg "Done: ${sample_id}"
}
export -f trim_one_sample

# ---- Run parallel trimming (2 concurrent × 2 threads = 4 cores) ------------
log_msg "Trimming samples (${PARALLEL_JOBS} in parallel, 2 threads each)..."

tail -n +2 "${SAMPLES_TSV}" | cut -f1 | \
    parallel -j "${PARALLEL_JOBS}" --halt soon,fail=1 \
    trim_one_sample {}

# ---- Summary report ---------------------------------------------------------
log_msg "=== fastp QC summary ==="
echo -e "sample\ttotal_reads\tfiltered_reads\tq20_rate\tq30_rate\tadapter_rate" > "${QC_DIR}/qc_summary.tsv"

for json_file in "${QC_DIR}"/*_fastp.json; do
    [[ -f "${json_file}" ]] || continue
    sample=$(basename "${json_file}" _fastp.json)
    python3 -c "
import json, sys
with open('${json_file}') as f:
    d = json.load(f)
s = d['summary']
bf = s['before_filtering']
af = s['after_filtering']
adapt = d.get('adapter_cutting', {})
adapt_rate = adapt.get('adapter_trimmed_reads', 0) / max(bf['total_reads'], 1) * 100
# q20_rate and q30_rate in fastp are floats 0-1 (proportion, not percentage)
q20 = af['q20_rate']
q30 = af['q30_rate']
# Handle both formats: some fastp versions output as ratio, others as percentage
if q20 <= 1:
    q20 *= 100
    q30 *= 100
print(f'${sample}\t{bf[\"total_reads\"]}\t{af[\"total_reads\"]}\t{q20:.1f}%\t{q30:.1f}%\t{adapt_rate:.1f}%')
" >> "${QC_DIR}/qc_summary.tsv"
done

log_msg "QC summary saved to ${QC_DIR}/qc_summary.tsv"
cat "${QC_DIR}/qc_summary.tsv" | column -t
log_msg "=== Step 3 complete ==="
