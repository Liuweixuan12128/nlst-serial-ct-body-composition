#!/bin/bash
#SBATCH --job-name=p2maps
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/p2_maps.%A_%a.out
P2_PROJECT_ROOT="${P2_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$P2_PROJECT_ROOT"
mkdir -p logs
./env/bin/python pipeline/p2_maps.py --out results/p2/maps --shard ${SLURM_ARRAY_TASK_ID}/4 --workers 44
