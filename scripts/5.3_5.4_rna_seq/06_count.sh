#!/usr/bin/env bash
# =============================================================================
# 06_count.sh — Feature counting with featureCounts + matrix generation
#
# FIXES applied:
#   - CRITICAL: heredoc delimiter unquoted for bash variable interpolation
#   - Added -B (require both mates), --primary, -Q 10 to featureCounts
#   - Dynamic outlier exclusion based on OUTLIER_HANDLING config
#   - Uses process substitution instead of pipe|while
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

log_msg "=== Step 6: Feature counting with featureCounts ==="

# ---- Determine strandedness -------------------------------------------------
if [[ "${STRAND}" == "auto" ]]; then
    DETECTED_FILE="${SALMON_CHECK_DIR}/detected_strand.txt"
    if [[ -f "${DETECTED_FILE}" ]]; then
        STRAND=$(cat "${DETECTED_FILE}")
        log_msg "Using auto-detected strandedness: -s ${STRAND}"
    else
        log_msg "WARNING: No strandedness detection file found. Using unstranded (-s 0)."
        STRAND="0"
    fi
fi

# ---- Verify GTF exists ------------------------------------------------------
if [[ ! -f "${GTF}" ]]; then
    log_msg "ERROR: GTF file not found at ${GTF}. Run 01_prepare_refs.sh first."
    exit 1
fi

# ---- Collect BAM files in sample order (with outlier handling) ---------------
> "${COUNTS_DIR}/bam_list.txt"
> "${COUNTS_DIR}/sample_order.txt"

while IFS=$'\t' read -r sample_id srr_id tissue condition replicate outlier_flag; do
    # Dynamic outlier exclusion
    if [[ "${OUTLIER_HANDLING}" == "exclude" && "${sample_id}" == "${OUTLIER_SAMPLE}" ]]; then
        log_msg "EXCLUDING outlier sample: ${sample_id} (${OUTLIER_SRR})"
        continue
    fi
    echo "${BAM_DIR}/${sample_id}.bam" >> "${COUNTS_DIR}/bam_list.txt"
    echo "${sample_id}" >> "${COUNTS_DIR}/sample_order.txt"
done < <(tail -n +2 "${SAMPLES_TSV}")

# Verify all BAM files exist
while read -r bam_path; do
    if [[ ! -f "${bam_path}" ]]; then
        log_msg "ERROR: BAM file not found: ${bam_path}"
        exit 1
    fi
done < "${COUNTS_DIR}/bam_list.txt"

BAM_FILES=$(cat "${COUNTS_DIR}/bam_list.txt" | tr '\n' ' ')

# ---- Run featureCounts -------------------------------------------------------
log_msg "Running featureCounts on all samples..."
featureCounts \
    -F GTF \
    -t "${FC_FEATURE_TYPE}" \
    -g "${FC_ATTR_TYPE}" \
    -p --countReadPairs \
    -B \
    --primary \
    -Q 10 \
    -s "${STRAND}" \
    -T "${THREADS}" \
    -a "${GTF}" \
    -o "${COUNTS_DIR}/all_counts.txt" \
    ${BAM_FILES} \
    2>&1 | tee "${LOG_DIR}/featureCounts.log"

log_msg "featureCounts complete."

# ---- Parse featureCounts output into clean CSV matrices ----------------------
log_msg "Parsing featureCounts output into clean CSV matrices..."

# FIX: No quotes on heredoc delimiter — bash variables MUST be interpolated
python3 << PYEOF
import pandas as pd
import os
import sys

counts_dir = "${COUNTS_DIR}"
design_dir = "${DESIGN_DIR}"

# Read featureCounts output (skip first comment line)
fc = pd.read_csv(f"{counts_dir}/all_counts.txt", sep="\t", comment="#")

# The columns are: Geneid, Chr, Start, End, Strand, Length, then BAM paths
# Rename BAM path columns to sample IDs
bam_cols = [c for c in fc.columns if c.endswith(".bam")]
sample_names = []
with open(f"{counts_dir}/sample_order.txt") as f:
    sample_names = [line.strip() for line in f if line.strip()]

if len(bam_cols) != len(sample_names):
    print(f"ERROR: {len(bam_cols)} BAM columns but {len(sample_names)} sample names")
    sys.exit(1)

rename_map = dict(zip(bam_cols, sample_names))
fc = fc.rename(columns=rename_map)

# Save gene lengths for TPM calculation
gene_lengths = fc[["Geneid", "Length"]].copy()
gene_lengths.columns = ["gene_id", "length"]
gene_lengths.to_csv(f"{counts_dir}/gene_lengths.tsv", sep="\t", index=False)

# Extract count matrix (Geneid + sample columns)
count_cols = ["Geneid"] + sample_names
counts_all = fc[count_cols].copy()
counts_all = counts_all.rename(columns={"Geneid": "gene_id"})
counts_all = counts_all.set_index("gene_id")
counts_all.to_csv(f"{counts_dir}/counts_all.csv")
print(f"Full count matrix: {counts_all.shape[0]} genes x {counts_all.shape[1]} samples")

# --- Generate subset matrices for each analysis ---

# DE Leaf: LeafCtrl + LeafStress samples
leaf_samples = [s for s in sample_names if s.startswith("Leaf")]
counts_leaf = counts_all[leaf_samples]
counts_leaf.to_csv(f"{counts_dir}/counts_leaf_de.csv")
print(f"Leaf DE matrix: {counts_leaf.shape}")

# DE Roots: RootsCtrl + RootsStress samples
roots_samples = [s for s in sample_names if s.startswith("Roots")]
counts_roots = counts_all[roots_samples]
counts_roots.to_csv(f"{counts_dir}/counts_roots_de.csv")
print(f"Roots DE matrix: {counts_roots.shape}")

# Basal: all control/basal samples
basal_design = pd.read_csv(f"{design_dir}/design_basal.csv")
basal_samples = basal_design["sample"].tolist()
counts_basal = counts_all[[s for s in basal_samples if s in counts_all.columns]]
counts_basal.to_csv(f"{counts_dir}/counts_basal.csv")
print(f"Basal matrix: {counts_basal.shape}")

print("Done parsing count matrices.")
PYEOF

# ---- Verify aquaporin genes are in count matrix ------------------------------
log_msg "Checking aquaporin gene coverage..."

# FIX: No quotes on heredoc delimiter
python3 << PYEOF
import pandas as pd

counts = pd.read_csv("${COUNTS_DIR}/counts_all.csv", index_col=0)
aqp = pd.read_csv("${AQUAPORIN_TABLE}", sep="\t")
aqp_ids = set(aqp["gene_id"].tolist())
found = aqp_ids.intersection(set(counts.index))
missing = aqp_ids - found

print(f"Aquaporin genes in count matrix: {len(found)}/{len(aqp_ids)}")
if missing:
    print(f"WARNING: Missing aquaporin genes: {missing}")
else:
    print("All aquaporin genes found in count matrix.")
PYEOF

log_msg "=== Step 6 complete ==="
