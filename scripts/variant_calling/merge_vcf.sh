#!/usr/bin/env bash
#SBATCH --job-name=MergeVCFs
#SBATCH --output=logs/MergeVCFs_%j.out
#SBATCH --error=logs/MergeVCFs_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00

set -euo pipefail

CONFIG="$PWD/config.yaml"
SAMPLESHEET="$PWD/samplesheet.tsv"

get_yaml() { python3 "$PWD/get_yaml.py" "$CONFIG" "$1"; }

FILTERED_DIR=$(get_yaml filtered_dir)
MERGED_DIR="$PWD/results/Merged"
mkdir -p "$MERGED_DIR" logs

echo "[$(date)] MergeVCFs - starting merge per sample..." | tee -a logs/MergeVCFs_master.log

# For each sample
tail -n +2 "$SAMPLESHEET" | cut -f1 | sort -u | while read -r SAMPLE; do
    echo ">> Processing sample: $SAMPLE" | tee -a logs/MergeVCFs_master.log

    mapfile -t files < <(find "$FILTERED_DIR" -maxdepth 1 -type f -name "${SAMPLE}_*_filtered.vcf.gz" -print | sort -V)
    OUTPUT_MERGED="${MERGED_DIR}/${SAMPLE}_merged.vcf.gz"
    LIST_PATH="${MERGED_DIR}/${SAMPLE}_vcfs.list"

    # Idempotency
    if [[ -s "$OUTPUT_MERGED" ]]; then
        echo "  - Merge already exists: $OUTPUT_MERGED — skipping" | tee -a logs/MergeVCFs_master.log
        continue
    fi

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "  - No VCF found for $SAMPLE — skipping" | tee -a logs/MergeVCFs_master.log
        continue
    fi

    printf "%s\n" "${files[@]}" > "$LIST_PATH"
    echo ">> Input list: $LIST_PATH" | tee -a logs/MergeVCFs_master.log

    # Execute merge
    if command -v picard >/dev/null 2>&1; then
        picard MergeVcfs I="@${LIST_PATH}" O="$OUTPUT_MERGED" &> "logs/MergeVCFs_${SAMPLE}.log"
        RET=$?
    else
        echo "ERROR: picard not found in PATH; adjust the script" | tee -a logs/MergeVCFs_master.log >&2
        RET=10
    fi

    if [[ $RET -ne 0 ]]; then
        echo "ERROR merging VCFs for $SAMPLE — see logs/MergeVCFs_${SAMPLE}.log" | tee -a logs/MergeVCFs_master.log >&2
        continue
    fi

    echo "Merge completed for $SAMPLE: $OUTPUT_MERGED" | tee -a logs/MergeVCFs_master.log
done

echo "[$(date)] MergeVCFs finished for all samples." | tee -a logs/MergeVCFs_master.log
