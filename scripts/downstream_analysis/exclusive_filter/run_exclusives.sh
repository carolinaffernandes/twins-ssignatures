#!/bin/bash
#SBATCH --job-name=exclusive_annovar
#SBATCH --output=logs/exclusive_annovar_%j.out
#SBATCH --error=logs/exclusive_annovar_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# conda activate variant-calling

echo "Initializing pipeline at $(date)"

python3 exclusives.py

echo "Initializing finished at $(date)"

