#!/usr/bin/env bash
# =============================================================================
# 02b_download_galaxy_bams.sh — Download pre-aligned BAMs from Galaxy
#
# BAMs are stored on HDD for space efficiency.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
check_env
ensure_dirs

log_msg "=== Step 2B: Downloading pre-aligned BAMs from Galaxy ==="

# Galaxy BAM URLs mapped to sample IDs
declare -A GALAXY_URLS=(
    ["GreenFruit_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135829e7652cfa5ec84/display?to_ext=bam"
    ["GreenFruit_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a47213589823cee72501f75/display?to_ext=bam"
    ["GreenFruit_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135ed4e23d50511f115/display?to_ext=bam"
    ["RedFruit_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a47213564dd9ca53a39a871/display?to_ext=bam"
    ["RedFruit_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a4721354cb2af9e4e5e89a6/display?to_ext=bam"
    ["RedFruit_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135d2fcab341bd3a4bc/display?to_ext=bam"
    ["Crown_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a4721355f3a161e6caaf9e7/display?to_ext=bam"
    ["Crown_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a47213508617a2c5a5da0a5/display?to_ext=bam"
    ["Crown_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a4721355fdffcb22e8ad02f/display?to_ext=bam"
    ["LeafStress_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135daa23980e2b796b7/display?to_ext=bam"
    ["LeafStress_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135d044406931bb99f1/display?to_ext=bam"
    ["LeafStress_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135b70f193bdf00f0fb/display?to_ext=bam"
    ["LeafCtrl_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135b1677a8362e0eec4/display?to_ext=bam"
    ["LeafCtrl_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a47213565df2122589f9017/display?to_ext=bam"
    ["LeafCtrl_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135028b757242394c31/display?to_ext=bam"
    ["RootsStress_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135cefa4bac48b22e02/display?to_ext=bam"
    ["RootsStress_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a47213560517cc22d901077/display?to_ext=bam"
    ["RootsStress_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a4721352d00254ce8f63f37/display?to_ext=bam"
    ["RootsCtrl_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135092ce7cf39c0cab7/display?to_ext=bam"
    ["RootsCtrl_2"]="https://usegalaxy.org/api/datasets/f9cad7b01a4721352f4cfb0b50fe9c48/display?to_ext=bam"
    ["RootsCtrl_3"]="https://usegalaxy.org/api/datasets/f9cad7b01a472135d41c3be683f190dc/display?to_ext=bam"
    ["AuxBud_1"]="https://usegalaxy.org/api/datasets/f9cad7b01a4721359d8db19290506906/display?to_ext=bam"
)

for sample_id in "${!GALAXY_URLS[@]}"; do
    BAM_FILE="${BAM_DIR}/${sample_id}.bam"
    BAM_INDEX="${BAM_FILE}.bai"

    if [[ -f "${BAM_FILE}" && -f "${BAM_INDEX}" ]]; then
        log_msg "Skipping ${sample_id} — already downloaded."
        continue
    fi

    log_msg "Downloading ${sample_id} from Galaxy..."
    wget -q --show-progress -O "${BAM_FILE}" "${GALAXY_URLS[${sample_id}]}" \
        2>&1 | tee -a "${LOG_DIR}/${sample_id}_galaxy_download.log" || {
        log_msg "ERROR: Failed to download ${sample_id}. Galaxy URL may have expired."
        rm -f "${BAM_FILE}"
        continue
    }

    log_msg "Indexing ${sample_id}..."
    samtools index "${BAM_FILE}"
    log_msg "Done: ${sample_id}"
done

log_msg "=== Step 2B complete ==="
log_msg "Downloaded BAM files:"
ls -lh "${BAM_DIR}"/*.bam 2>/dev/null || log_msg "No BAM files found!"
