#!/usr/bin/env bash
#SBATCH --job-name=run_pipeline
#SBATCH --output=logs/run_pipeline_%j.out
#SBATCH --error=logs/run_pipeline_%j.err
#SBATCH --time=192:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

set -euo pipefail

# ===============================================
# PARÂMETROS
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

get_yaml() { python3 "$PWD/get_yaml.py" "$CONFIG" "$1"; }

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

echo "=== CONFIGURAÇÕES CARREGADAS ==="
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
# GERAR TASK LIST (amostra × cromossomo)
# ===============================================
echo ">> Gerando lista de tarefas (amostra × cromossomo)..."
> "$TASK_LIST"
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE BAM_PATH; do
    while read -r CHROM; do
        echo -e "${SAMPLE}\t${BAM_PATH}\t${CHROM}" >> "$TASK_LIST"
    done < "$CHROM_FILE"
done

TOTAL_TASKS=$(wc -l < "$TASK_LIST")
echo "Arquivo $TASK_LIST criado com $TOTAL_TASKS combinações."

# ===============================================
# FUNÇÃO PARA SUBMISSÃO EM CHUNKS
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

        echo "[DEBUG] Submetido bloco ${start}-${end} ($step) JOBID=$JOBID"
        jobids+=("$JOBID")

        # Aguarda o bloco terminar
        echo "[DEBUG] Aguardando bloco ${start}-${end} de ${step}..."
        while squeue -j "$JOBID" 2>/dev/null | grep -q "$USER"; do
            sleep 30
            if [[ $DEBUG -eq 1 ]]; then
                echo "[DEBUG] Ainda rodando tasks do bloco ${start}-${end}..."
                squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"
            fi
        done

        start=$((end + 1))
    done

    # Retorna todos os jobids concatenados por vírgula para dependência correta
    echo "$(IFS=,; echo "${jobids[*]}")"
}

# ===============================================
# SUBMISSÃO
# ===============================================
echo "[DEBUG] Submetendo Mutect2..."
JOBID_MUTECT2=$(submit_in_chunks "scripts/mutect2_chr.sh" "$TOTAL_TASKS" "Mutect2" "")

echo "[DEBUG] Submetendo FilterMutectCalls..."
JOBID_FMC=$(submit_in_chunks "scripts/filtermutectcalls_chr.sh" "$TOTAL_TASKS" "FilterMutectCalls" "$JOBID_MUTECT2")

echo "[DEBUG] Submetendo MergeVCFs..."
JOBID_MERGE=$(sbatch --parsable --dependency=afterok:${JOBID_FMC} scripts/merge_vcfs.sh)

echo "=============================================="
echo "Pipeline submetida com sucesso!"
echo "Mutect2 JOBID final: $JOBID_MUTECT2"
echo "FilterMutectCalls JOBID final: $JOBID_FMC"
echo "MergeVCFs JOBID: $JOBID_MERGE"
echo "=============================================="

