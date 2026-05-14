#!/bin/bash
#SBATCH --job-name=exclusive_annovar
#SBATCH --output=logs/exclusive_annovar_%j.out
#SBATCH --error=logs/exclusive_annovar_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# OPCIONAL: ativar ambiente Conda ou módulo Python
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate variant-calling

echo "Iniciando pipeline às $(date)"

python3 exclusivas.py

echo "Pipeline finalizado às $(date)"

