#!/usr/bin/env bash

REF_DIR="$PWD/ref"
mkdir -p "$REF_DIR"

echo ">> Baixando referências para $REF_DIR ..."

# Referência do genoma hg38
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta* "$REF_DIR/"

# Painel de normais (PON)
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz* "$REF_DIR/"

# Germline resource
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz* "$REF_DIR/"

echo ">> Downloads concluídos com sucesso!"
echo "Arquivos salvos em: $REF_DIR"

