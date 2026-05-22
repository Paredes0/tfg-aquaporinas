#!/usr/bin/env bash
# =============================================================================
# 05_strandedness_check.sh — Auto-detect library strandedness using salmon
#
# FIXES applied:
#   - CRITICAL: heredoc delimiter unquoted so bash variables are interpolated
#   - Uses process substitution instead of pipe|while
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

log_msg "=== Step 5: Strandedness detection ==="

if [[ "${STRAND}" != "auto" ]]; then
    log_msg "Strandedness manually set to: ${STRAND}. Skipping auto-detection."
    exit 0
fi

# Verify salmon index exists
if [[ ! -d "${SALMON_IDX}" ]]; then
    log_msg "ERROR: Salmon index not found. Run 01_prepare_refs.sh first."
    exit 1
fi

# Pick one representative sample per experiment to test strandedness
declare -a TEST_SAMPLES=("GreenFruit_1" "LeafCtrl_1" "RootsCtrl_1" "Crown_1" "AuxBud_1")

echo -e "sample\tlib_type\tcompatible_frags\tISF\tISR\tIU" > "${SALMON_CHECK_DIR}/strandedness_results.tsv"

for sample_id in "${TEST_SAMPLES[@]}"; do

    R1="${TRIMMED_DIR}/${sample_id}_1.fastq.gz"
    R2="${TRIMMED_DIR}/${sample_id}_2.fastq.gz"
    SALMON_OUT="${SALMON_CHECK_DIR}/${sample_id}"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        # If using Galaxy BAMs, use a subset of reads from the BAM
        BAM_FILE="${BAM_DIR}/${sample_id}.bam"
        if [[ -f "${BAM_FILE}" ]]; then
            log_msg "Extracting subset from BAM for ${sample_id}..."
            mkdir -p "${SALMON_CHECK_DIR}/temp_fq"
            samtools fastq -@ "${SORT_THREADS}" \
                -1 "${SALMON_CHECK_DIR}/temp_fq/${sample_id}_1.fq.gz" \
                -2 "${SALMON_CHECK_DIR}/temp_fq/${sample_id}_2.fq.gz" \
                -0 /dev/null -s /dev/null \
                --reference "${GENOME}" \
                "${BAM_FILE}" 2>/dev/null || true
            R1="${SALMON_CHECK_DIR}/temp_fq/${sample_id}_1.fq.gz"
            R2="${SALMON_CHECK_DIR}/temp_fq/${sample_id}_2.fq.gz"
        else
            log_msg "WARNING: No reads available for ${sample_id}, skipping."
            continue
        fi
    fi

    if [[ -d "${SALMON_OUT}" ]]; then
        log_msg "Skipping ${sample_id} — already checked."
    else
        log_msg "Running salmon on ${sample_id} (auto-detect mode)..."
        salmon quant \
            -i "${SALMON_IDX}" \
            -l A \
            -1 "${R1}" -2 "${R2}" \
            --validateMappings \
            --threads "${THREADS}" \
            -o "${SALMON_OUT}" \
            2>&1 | tee "${LOG_DIR}/${sample_id}_salmon_check.log"
    fi

    # Parse the library format counts
    LIB_JSON="${SALMON_OUT}/lib_format_counts.json"
    if [[ -f "${LIB_JSON}" ]]; then
        # FIX: No quotes on heredoc delimiter — bash variables MUST be interpolated
        python3 << PYEOF
import json
with open('${LIB_JSON}') as f:
    d = json.load(f)
expected = d.get('expected_format', 'unknown')
compat = d.get('compatible_fragment_ratio', 0)
# In salmon 1.10.3, ISF/ISR/IU are top-level keys (not nested under read_files)
isf = d.get('ISF', 0)
isr = d.get('ISR', 0)
iu  = d.get('IU', 0)
print(f'${sample_id}\t{expected}\t{compat:.4f}\t{isf}\t{isr}\t{iu}')
PYEOF
    fi
done >> "${SALMON_CHECK_DIR}/strandedness_results.tsv"

# ---- Determine consensus strandedness ---------------------------------------
log_msg "=== Strandedness results ==="
cat "${SALMON_CHECK_DIR}/strandedness_results.tsv" | column -t

# FIX: No quotes on heredoc delimiter — bash variables MUST be interpolated
python3 << PYEOF
import sys

results_file = "${SALMON_CHECK_DIR}/strandedness_results.tsv"
strand_map = {"ISF": "1", "ISR": "2", "IU": "0"}
detected = []

with open(results_file) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            lib_type = parts[1]
            detected.append(lib_type)

if not detected:
    print("WARNING: No strandedness detected. Defaulting to unstranded (0).")
    strand = "0"
else:
    # Map salmon lib types to featureCounts -s values
    strand_vals = []
    for lt in detected:
        if "SR" in lt or lt == "ISR":
            strand_vals.append("2")
        elif "SF" in lt or lt == "ISF":
            strand_vals.append("1")
        else:
            strand_vals.append("0")

    # Use majority vote
    from collections import Counter
    consensus = Counter(strand_vals).most_common(1)[0][0]
    strand = consensus

    label = {"0": "unstranded", "1": "stranded (sense)", "2": "reversely stranded"}
    print(f"Consensus strandedness: {label.get(strand, strand)} (featureCounts -s {strand})")

# Write the detected strandedness to a file for downstream scripts
with open("${SALMON_CHECK_DIR}/detected_strand.txt", "w") as f:
    f.write(strand)
PYEOF

log_msg "Strandedness value saved to ${SALMON_CHECK_DIR}/detected_strand.txt"
log_msg "=== Step 5 complete ==="
