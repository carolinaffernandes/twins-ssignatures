#!/usr/bin/env bash
#SBATCH --job-name=FilterMutectCalls
#SBATCH --output=logs/FilterMutectCalls_%A_%a.out
#SBATCH --error=logs/FilterMutectCalls_%A_%a.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00

set -euo pipefail

CONFIG="$PWD/config.yaml"
TASK_LIST="$PWD/task_list.txt"

DEBUG=0
if [[ "${1:-}" == "--debug" ]]; then
    DEBUG=1
    set -x
    echo "=== DEBUG MODE ENABLED ==="
fi

get_yaml() { python3 "$PWD/get_yaml.py" "$CONFIG" "$1"; }

OUTPUT_DIR=$(get_yaml output_dir)
FILTERED_DIR=$(get_yaml filtered_dir)
REFERENCE=$(get_yaml reference)
SINGULARITY_IMG=$(get_yaml singularity_img)

TASK_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")
SAMPLE=$(echo "$TASK_LINE" | cut -f1)
CHROM=$(echo "$TASK_LINE" | cut -f3)

MUTECT2_OUT="${OUTPUT_DIR}/${SAMPLE}_${CHROM}.vcf.gz"
FILTERED_OUT="${FILTERED_DIR}/${SAMPLE}_${CHROM}_filtered.vcf.gz"
LOG_DIR="logs"
LOG_FILE="${LOG_DIR}/${SAMPLE}_FilterMutectCalls_${CHROM}.log"

mkdir -p "$FILTERED_DIR" "$LOG_DIR"

echo "[$(date)] Running FilterMutectCalls for $SAMPLE $CHROM" | tee -a "$LOG_FILE"
echo "Input: $MUTECT2_OUT" | tee -a "$LOG_FILE"
echo "Output: $FILTERED_OUT" | tee -a "$LOG_FILE"

# === Idempotência ===
if [[ -s "$FILTERED_OUT" ]]; then
    echo "[$(date)] SKIP: $FILTERED_OUT already exists" | tee -a "$LOG_FILE"
    exit 0
fi

if [[ ! -s "$MUTECT2_OUT" ]]; then
    echo "[$(date)] ERROR: input VCF missing $MUTECT2_OUT" | tee -a "$LOG_FILE" >&2
    exit 2
fi

# Executa
CMD="singularity exec $SINGULARITY_IMG gatk FilterMutectCalls -R $REFERENCE -V $MUTECT2_OUT -O $FILTERED_OUT"
echo "Command: $CMD" | tee -a "$LOG_FILE"
eval "$CMD" &>> "$LOG_FILE"

if [[ $? -ne 0 ]]; then
    echo "[$(date)] ERROR in FilterMutectCalls for $SAMPLE $CHROM" | tee -a "$LOG_FILE" >&2
    exit 1
fi

echo "DONE FILTERMUTECTCALLS" >> "$LOG_FILE"
echo "[$(date)] FilterMutectCalls finished for $SAMPLE $CHROM" | tee -a "$LOG_FILE"

