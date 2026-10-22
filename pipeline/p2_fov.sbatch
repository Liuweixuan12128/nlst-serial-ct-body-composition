#!/bin/bash
#SBATCH --job-name=p2fov
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --output=logs/p2_fov.%j.out
P2_PROJECT_ROOT="${P2_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$P2_PROJECT_ROOT"
mkdir -p logs
./env/bin/python pipeline/p2_fov_truncation.py --out results/p2/p2_fov_truncation.csv --workers 44
