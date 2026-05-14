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

echo "[$(date)] MergeVCFs - iniciando merge por amostra..." | tee -a logs/MergeVCFs_master.log

# Para cada amostra
tail -n +2 "$SAMPLESHEET" | cut -f1 | sort -u | while read -r SAMPLE; do
    echo ">> Processando amostra: $SAMPLE" | tee -a logs/MergeVCFs_master.log

    mapfile -t files < <(find "$FILTERED_DIR" -maxdepth 1 -type f -name "${SAMPLE}_*_filtered.vcf.gz" -print | sort -V)
    OUTPUT_MERGED="${MERGED_DIR}/${SAMPLE}_merged.vcf.gz"
    LIST_PATH="${MERGED_DIR}/${SAMPLE}_vcfs.list"

    # Idempotência
    if [[ -s "$OUTPUT_MERGED" ]]; then
        echo "  - Merge já existe: $OUTPUT_MERGED — pulando" | tee -a logs/MergeVCFs_master.log
        continue
    fi

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "  - Nenhum VCF encontrado para $SAMPLE — pulando" | tee -a logs/MergeVCFs_master.log
        continue
    fi

    printf "%s\n" "${files[@]}" > "$LIST_PATH"
    echo ">> Lista de entrada: $LIST_PATH" | tee -a logs/MergeVCFs_master.log

    # Executa merge
    if command -v picard >/dev/null 2>&1; then
        picard MergeVcfs I="@${LIST_PATH}" O="$OUTPUT_MERGED" &> "logs/MergeVCFs_${SAMPLE}.log"
        RET=$?
    else
        echo "ERRO: picard não encontrado no PATH; ajuste o script" | tee -a logs/MergeVCFs_master.log >&2
        RET=10
    fi

    if [[ $RET -ne 0 ]]; then
        echo "ERRO ao unir VCFs de $SAMPLE — veja logs/MergeVCFs_${SAMPLE}.log" | tee -a logs/MergeVCFs_master.log >&2
        continue
    fi

    echo "Merge concluído para $SAMPLE: $OUTPUT_MERGED" | tee -a logs/MergeVCFs_master.log
done

echo "[$(date)] MergeVCFs finalizado para todas as amostras." | tee -a logs/MergeVCFs_master.log

