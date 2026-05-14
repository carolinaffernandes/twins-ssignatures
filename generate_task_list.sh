#!/bin/bash
SAMPLESHEET="$PWD/samplesheet.tsv"
CHROM_FILE="$PWD/chrom_list.txt"
TASK_LIST="$PWD/task_list.txt"

> "$TASK_LIST"
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r SAMPLE BAM_PATH; do
  for CHROM in $(cat "$CHROM_FILE"); do
    echo -e "${SAMPLE}\t${BAM_PATH}\t${CHROM}" >> "$TASK_LIST"
  done
done

echo "Arquivo $TASK_LIST gerado com $(wc -l < "$TASK_LIST") combinações."

