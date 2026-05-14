#!/usr/bin/env bash
#SBATCH --job-name=Mutect2
#SBATCH --output=logs/Mutect2_%A_%a.out
#SBATCH --error=logs/Mutect2_%A_%a.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00

set -euo pipefail

CONFIG="$PWD/config.yaml"
TASK_LIST="$PWD/task_list.txt"

get_yaml() { python3 get_yaml.py "$CONFIG" "$1"; }

OUTPUT_DIR=$(get_yaml output_dir)
REFERENCE=$(get_yaml reference)
GERMLINE=$(get_yaml germline_resource)
PON=$(get_yaml pon)
SINGULARITY_IMG=$(get_yaml singularity_img)

TASK_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")
SAMPLE=$(echo "$TASK_LINE" | cut -f1)
BAM_PATH=$(echo "$TASK_LINE" | cut -f2)
CHROM=$(echo "$TASK_LINE" | cut -f3)

MUTECT2_OUT="${OUTPUT_DIR}/${SAMPLE}_${CHROM}.vcf.gz"
LOG_DIR="logs"
LOG_FILE="${LOG_DIR}/${SAMPLE}_Mutect2_${CHROM}.log"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Se já existe VCF final não reprocessa
if [[ -s "$MUTECT2_OUT" ]]; then
  echo "[$(date)] SKIP: ${MUTECT2_OUT} already exists — skipping Mutect2." | tee -a "$LOG_FILE"
  exit 0
fi

# Se log contém marca de conclusão, pula
if grep -q "^DONE MUTECT2$" "$LOG_FILE" 2>/dev/null; then
  echo "[$(date)] SKIP: ${LOG_FILE} contains DONE MUTECT2 — skipping." | tee -a "$LOG_FILE"
  exit 0
fi

echo ">> $SAMPLE - Mutect2 - $CHROM" | tee -a "$LOG_FILE"
singularity exec "$SINGULARITY_IMG" gatk Mutect2 \
    -R "$REFERENCE" \
    -I "$BAM_PATH" \
    --germline-resource "$GERMLINE" \
    --panel-of-normals "$PON" \
    -L "$CHROM" \
    -O "$MUTECT2_OUT" \
    &>> "$LOG_FILE"

if [[ $? -ne 0 ]]; then
    echo "[$(date)] ERRO no Mutect2 para $SAMPLE - $CHROM" | tee -a "$LOG_FILE" >&2
    exit 1
fi

echo "DONE MUTECT2" >> "$LOG_FILE"
echo "[$(date)] Mutect2 finished successfully for ${SAMPLE}_${CHROM}" | tee -a "$LOG_FILE"
