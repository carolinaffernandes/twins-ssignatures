#!/usr/bin/env bash
#SBATCH --job-name=run_pipeline
#SBATCH --output=logs/run_pipeline_%j.out
#SBATCH --error=logs/run_pipeline_%j.err
#SBATCH --time=192:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

set -euo pipefail

# ===============================================
# PARAMETERS
# ===============================================
CONFIG="$PWD/config.yaml"
SAMPLESHEET="$PWD/samplesheet.tsv"
CHROM_FILE="$PWD/chrom_list.txt"
TASK_LIST="$PWD/task_list.txt"


DEBUG=0
if [[ "${1:-}" == "--debug" ]]; then
    DEBUG=1
    set -x
    echo "=== DEBUG MODE ENABLED ==="
fi

get_yaml() { python3 "$PWD/scripts/utils/get_yaml.py" "$CONFIG" "$1"; }

MAX_SUBMIT=$(get_yaml max_submit)
MAX_JOBS=$(get_yaml max_jobs)
OUTPUT_DIR=$(get_yaml output_dir)
FILTERED_DIR=$(get_yaml filtered_dir)
TMP_DIR=$(get_yaml tmp_dir)
REFERENCE=$(get_yaml reference)
GERMLINE=$(get_yaml germline_resource)
PON=$(get_yaml pon)
SINGULARITY_IMG=$(get_yaml singularity_img)

mkdir -p "$OUTPUT_DIR" "$FILTERED_DIR" "$TMP_DIR" "logs"

echo "=== LOADED CONFIGURATIONS ==="
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "FILTERED_DIR: $FILTERED_DIR"
echo "TMP_DIR: $TMP_DIR"
echo "REFERENCE: $REFERENCE"
echo "GERMLINE: $GERMLINE"
echo "PON: $PON"
echo "SINGULARITY_IMG: $SINGULARITY_IMG"
echo "MAX_JOBS: $MAX_JOBS"
echo "MAX_SUBMIT (QOS limit): $MAX_SUBMIT"
echo "==============================="

# ===============================================
# GENERATE TASK LIST (sample × chromosome)
# ===============================================
echo ">> Generating task list (sample × chromosome)..."
> "$TASK_LIST"
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE BAM_PATH; do
    while read -r CHROM; do
        echo -e "${SAMPLE}\t${BAM_PATH}\t${CHROM}" >> "$TASK_LIST"
    done < "$CHROM_FILE"
done

TOTAL_TASKS=$(wc -l < "$TASK_LIST")
echo "File $TASK_LIST created with $TOTAL_TASKS combinations."

# ===============================================
# FUNCTION FOR CHUNK SUBMISSION
# ===============================================
submit_in_chunks() {
    local script="$1"
    local total="$2"
    local step="$3"
    local dep="$4"

    local jobids=()
    local start=1

    while [[ $start -le $total ]]; do
        local end=$((start + MAX_JOBS - 1))
        [[ $end -gt $total ]] && end=$total

        if [[ -z "$dep" ]]; then
            JOBID=$(sbatch --parsable --array=${start}-${end}%${MAX_SUBMIT} "$script")
        else
            JOBID=$(sbatch --parsable --dependency=afterok:${dep} --array=${start}-${end}%${MAX_SUBMIT} "$script")
        fi

        echo "[DEBUG] Submitted block ${start}-${end} ($step) JOBID=$JOBID"
        jobids+=("$JOBID")

        # Wait for the block to finish
        echo "[DEBUG] Waiting for block ${start}-${end} of ${step}..."
        while squeue -j "$JOBID" 2>/dev/null | grep -q "$USER"; do
            sleep 30
            if [[ $DEBUG -eq 1 ]]; then
                echo "[DEBUG] Still running tasks from block ${start}-${end}..."
                squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
            fi
        done

        start=$((end + 1))
    done

    # Return all jobids concatenated by comma for correct dependency
    echo "$(IFS=,; echo "${jobids[*]}")"
}

# ===============================================
# SUBMISSION
# ===============================================
echo "[DEBUG] Submitting Mutect2..."
JOBID_MUTECT2=$(submit_in_chunks "scripts/variant_calling/mutect2_chr.sh" "$TOTAL_TASKS" "Mutect2" "")

echo "[DEBUG] Submitting FilterMutectCalls..."
JOBID_FMC=$(submit_in_chunks "scripts/variant_calling/filtermutectcalls_chr.sh" "$TOTAL_TASKS" "FilterMutectCalls" "$JOBID_MUTECT2")

echo "[DEBUG] Submitting MergeVCFs..."
JOBID_MERGE=$(sbatch --parsable --dependency=afterok:${JOBID_FMC} scripts/variant_calling/merge_vcf.sh)

echo "=============================================="
echo "Pipeline submitted successfully!"
echo "Final Mutect2 JOBID: $JOBID_MUTECT2"
echo "Final FilterMutectCalls JOBID: $JOBID_FMC"
echo "Final MergeVCFs JOBID: $JOBID_MERGE"
echo "=============================================="
