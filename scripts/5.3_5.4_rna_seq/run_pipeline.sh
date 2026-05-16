#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master script for Fragaria x ananassa RNA-seq pipeline
#
# Usage:
#   bash run_pipeline.sh              # Run full pipeline (SRA download)
#   bash run_pipeline.sh galaxy       # Use Galaxy BAMs instead of SRA
#   bash run_pipeline.sh from_step N  # Resume from step N (e.g., from_step 6)
#
# FIXES applied:
#   - Error propagation through pipes using PIPESTATUS
#   - Exports environment variables for R/Python scripts
#   - Added MultiQC integration (step 10)
#   - Added optional cleanup step
#   - Disk space check before starting
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ---- Parse arguments --------------------------------------------------------
DATA_MODE="${1:-sra}"      # "sra" or "galaxy"
START_STEP="${2:-0}"       # Step to start from (for resume)

if [[ "${DATA_MODE}" == "from_step" ]]; then
    START_STEP="${2:-0}"
    DATA_MODE="sra"
fi

# ---- Export environment variables for R/Python scripts -----------------------
export PROJECT_DIR
export OUTLIER_HANDLING
export OUTLIER_SAMPLE
export DATA_MODE

echo "=============================================================="
echo " Fragaria x ananassa Aquaporin RNA-seq Pipeline"
echo "=============================================================="
echo " Data mode:      ${DATA_MODE}"
echo " Start step:     ${START_STEP}"
echo " Environment:    ${ENV_NAME}"
echo " Project dir:    ${PROJECT_DIR}"
echo " SSD storage:    ${SSD_DIR}"
echo " HDD storage:    ${HDD_DIR}"
echo " Threads:        ${THREADS}"
echo " Parallel jobs:  ${PARALLEL_JOBS}"
echo " Outlier:        ${OUTLIER_HANDLING} (${OUTLIER_SAMPLE})"
echo "=============================================================="

# ---- Pre-flight checks ------------------------------------------------------
ensure_dirs

# Check disk space
SSD_FREE_GB=$(df --output=avail "${SSD_DIR}" | tail -1 | awk '{printf "%.0f", $1/1024/1024}')
HDD_FREE_GB=$(df --output=avail "${HDD_DIR}" | tail -1 | awk '{printf "%.0f", $1/1024/1024}')
RAM_GB=$(awk '/MemTotal/ {printf "%.0f", $2/1024/1024}' /proc/meminfo)

echo ""
echo " Hardware check:"
echo "   SSD free:  ${SSD_FREE_GB} GB (need ~20 GB for indices/results)"
echo "   HDD free:  ${HDD_FREE_GB} GB (need ~300 GB for FASTQ/BAM)"
echo "   RAM total: ${RAM_GB} GB (need ≥12 GB for HISAT2)"
echo "   CPU cores: $(nproc)"
echo ""

if [[ "${HDD_FREE_GB}" -lt 200 ]]; then
    echo "WARNING: HDD has <200 GB free. Pipeline may run out of space!"
    echo "Press Ctrl+C within 10 seconds to abort..."
    sleep 10
fi

if [[ "${RAM_GB}" -lt 12 ]]; then
    echo "WARNING: System has <12 GB RAM. HISAT2 index building may fail!"
    echo "Consider adding swap space before proceeding."
fi

# ---- Step runner with resume support + error propagation --------------------
run_step() {
    local step_num=$1
    local step_name=$2
    local step_script=$3

    if [[ ${step_num} -lt ${START_STEP} ]]; then
        echo "[SKIP] Step ${step_num}: ${step_name} (before start step)"
        return 0
    fi

    echo ""
    echo "=============================================================="
    echo "[STEP ${step_num}] ${step_name}"
    echo "=============================================================="

    # Run with proper error propagation
    bash "${SCRIPT_DIR}/${step_script}" 2>&1 | tee -a "${LOG_DIR}/pipeline.log"
    local exit_code=${PIPESTATUS[0]}

    if [[ ${exit_code} -ne 0 ]]; then
        echo "[FAILED] Step ${step_num}: ${step_name} (exit code ${exit_code})"
        echo "[$(date)] FAILED at step ${step_num}: ${step_name}" >> "${LOG_DIR}/pipeline.log"
        exit ${exit_code}
    fi

    echo "[DONE] Step ${step_num}: ${step_name}"
}

# ---- Run R scripts -----------------------------------------------------------
run_r_step() {
    local step_num=$1
    local step_name=$2
    local r_script=$3

    if [[ ${step_num} -lt ${START_STEP} ]]; then
        echo "[SKIP] Step ${step_num}: ${step_name} (before start step)"
        return 0
    fi

    echo ""
    echo "=============================================================="
    echo "[STEP ${step_num}] ${step_name}"
    echo "=============================================================="

    cd "${PROJECT_DIR}"
    Rscript "${SCRIPT_DIR}/${r_script}" 2>&1 | tee -a "${LOG_DIR}/pipeline.log"
    local exit_code=${PIPESTATUS[0]}

    if [[ ${exit_code} -ne 0 ]]; then
        echo "[FAILED] Step ${step_num}: ${step_name} (exit code ${exit_code})"
        echo "[$(date)] FAILED at step ${step_num}: ${step_name}" >> "${LOG_DIR}/pipeline.log"
        exit ${exit_code}
    fi

    echo "[DONE] Step ${step_num}: ${step_name}"
}

# ---- Run Python scripts ------------------------------------------------------
run_py_step() {
    local step_num=$1
    local step_name=$2
    local py_script=$3

    if [[ ${step_num} -lt ${START_STEP} ]]; then
        echo "[SKIP] Step ${step_num}: ${step_name} (before start step)"
        return 0
    fi

    echo ""
    echo "=============================================================="
    echo "[STEP ${step_num}] ${step_name}"
    echo "=============================================================="

    cd "${PROJECT_DIR}"
    python3 "${SCRIPT_DIR}/${py_script}" 2>&1 | tee -a "${LOG_DIR}/pipeline.log"
    local exit_code=${PIPESTATUS[0]}

    if [[ ${exit_code} -ne 0 ]]; then
        echo "[FAILED] Step ${step_num}: ${step_name} (exit code ${exit_code})"
        echo "[$(date)] FAILED at step ${step_num}: ${step_name}" >> "${LOG_DIR}/pipeline.log"
        exit ${exit_code}
    fi

    echo "[DONE] Step ${step_num}: ${step_name}"
}

# ==============================================================================
# PIPELINE EXECUTION
# ==============================================================================

# Initialize log
echo "[$(date)] Pipeline started — mode=${DATA_MODE}, start_step=${START_STEP}" >> "${LOG_DIR}/pipeline.log"

# Step 1: Prepare references (genome, GFF3→GTF, HISAT2 index, salmon index)
run_step 1 "Prepare references" "01_prepare_refs.sh"

# Step 2: Download data
if [[ "${DATA_MODE}" == "galaxy" ]]; then
    run_step 2 "Download Galaxy BAMs" "02b_download_galaxy_bams.sh"
else
    run_step 2 "Download FASTQ from SRA" "02a_download_sra.sh"

    # Step 3: QC and trimming (only for SRA path)
    run_step 3 "QC and trimming (fastp)" "03_qc_trim.sh"

    # Step 4: Alignment (only for SRA path)
    run_step 4 "HISAT2 alignment" "04_align.sh"
fi

# Step 5: Strandedness detection
run_step 5 "Strandedness detection" "05_strandedness_check.sh"

# Step 6: Feature counting
run_step 6 "Feature counting (featureCounts)" "06_count.sh"

# Step 7: Differential expression
run_r_step 7 "Differential expression (DESeq2)" "07_de_analysis.R"

# Step 8: Basal aquaporin expression
run_r_step 8 "Basal aquaporin expression" "08_basal_expression.R"

# Step 9: Reannotation flags
run_py_step 9 "Reannotation flag analysis" "09_reannotation_flags.py"

# Step 10: MultiQC report
if [[ 10 -ge ${START_STEP} ]]; then
    echo ""
    echo "=============================================================="
    echo "[STEP 10] MultiQC integrated report"
    echo "=============================================================="

    if command -v multiqc &> /dev/null; then
        multiqc \
            "${QC_DIR}" "${LOG_DIR}" "${COUNTS_DIR}" \
            -o "${RESULTS_DIR}/multiqc" \
            --force \
            --title "Fragaria x ananassa Aquaporin RNA-seq" \
            2>&1 | tee -a "${LOG_DIR}/pipeline.log"
        echo "[DONE] Step 10: MultiQC report at ${RESULTS_DIR}/multiqc/"
    else
        echo "[SKIP] Step 10: multiqc not installed. Run: micromamba install -c bioconda multiqc"
    fi
fi

# Step 11: Homeolog grouping (phylogenetic tree-based)
run_py_step 11 "Homeolog grouping from phylogenetic tree" "11_homeolog_grouping.py"

# Step 12: Homeolog expression analysis (3 levels: individual, collapsed, dominance)
run_r_step 12 "Homeolog expression analysis (3 levels)" "13_homeolog_expression.R"

# Step 13: Homeolog differential expression (DESeq2 grouped)
run_r_step 13 "Homeolog differential expression" "14_homeolog_de_analysis.R"

# Step 14 (OPTIONAL): GFF3 substitution for corrected aquaporin annotations
# This generates a new GFF3 with exonerate-corrected gene models.
# Re-run steps 1, 6-8, 11-13 using the corrected GFF3 to compare results.
# Uncomment to run:
# run_py_step 14 "GFF3 substitution (corrected AQP annotations)" "12_substitute_gff3.py"

# Step 15: Generate Homeolog eFP expression viewer
run_py_step 15 "Generate Homeolog eFP viewer" "15_homeolog_efp_viewer.py"

# ==============================================================================
# OPTIONAL CLEANUP — Uncomment to free disk space after successful pipeline run
# ==============================================================================
# echo ""
# echo "=============================================================="
# echo "[CLEANUP] Removing intermediate files to save disk space"
# echo "=============================================================="
# # Remove raw FASTQ (keep trimmed)
# rm -rf "${FASTQ_DIR}"/*.fastq.gz
# echo "Removed raw FASTQ files."
# # Remove SRA cache
# rm -rf "${SRA_CACHE_DIR}"
# echo "Removed SRA cache."
# # Remove temp files
# rm -rf "${TMP_DIR}"
# echo "Removed temp files."

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
echo ""
echo "=============================================================="
echo " PIPELINE COMPLETE"
echo "=============================================================="
echo ""
echo " Storage layout:"
echo "   SSD (fast):  ${SSD_DIR}"
echo "     - Indices, references, counts, results"
echo "   HDD (large): ${HDD_DIR}"
echo "     - FASTQ, trimmed reads, BAM files"
echo ""
echo " Results:"
echo "   DE Leaf:     ${RESULTS_DIR}/de_leaf/"
echo "   DE Roots:    ${RESULTS_DIR}/de_roots/"
echo "   Basal AQP:   ${RESULTS_DIR}/basal_aquaporins/"
echo "   Homeolog:    ${RESULTS_DIR}/homeolog_analysis/"
echo "   MultiQC:     ${RESULTS_DIR}/multiqc/"
echo ""
echo " Key output files:"
echo "   DE results:  results/de_leaf/results_leaf.csv"
echo "                results/de_roots/results_roots.csv"
echo "   AQP DE:      results/de_leaf/de_aquaporins_leaf.csv"
echo "                results/de_roots/de_aquaporins_roots.csv"
echo "   Basal TPM:   results/basal_aquaporins/basal_aquaporins_tpm.csv"
echo "   Reannotation: results/basal_aquaporins/reannotation_candidates.tsv"
echo "   IGV regions: results/basal_aquaporins/igv_aquaporin_regions.bed"
echo "   Homeolog:   homeolog_groups.tsv"
echo "   Collapsed:  results/homeolog_analysis/collapsed_tpm.csv"
echo "   Dominance:  results/homeolog_analysis/dominance_by_tissue.csv"
echo "   HG DE:      results/homeolog_de_analysis/de_homeologs_leaf.csv"
echo "   eFP Viewer: results/homeolog_analysis/efp_viewer_homeologs.html"
echo ""
echo " Plots:"
echo "   PCA, volcano, MA, heatmaps in results/de_leaf/ and results/de_roots/"
echo "   Basal heatmaps, profiles, barplots in results/basal_aquaporins/"
echo "   Homeolog: collapsed heatmap, dominance, PCA in results/homeolog_analysis/"
echo ""
echo " Notes:"
echo "   - AuxBud_1 (N=1) included as descriptive reference only"
echo "   - Outlier handling: ${OUTLIER_HANDLING} (${OUTLIER_SAMPLE})"
echo ""
echo " Logs: ${LOG_DIR}/"
echo "=============================================================="
echo "[$(date)] Pipeline completed" >> "${LOG_DIR}/pipeline.log"
