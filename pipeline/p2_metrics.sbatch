#!/bin/bash
#SBATCH --job-name=p2met
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --output=logs/p2_metrics.%j.out
# Paper 2: paraspinal measures for every segmented serial scan (resumes: scans already measured without error are kept).
P2_PROJECT_ROOT="${P2_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$P2_PROJECT_ROOT"
mkdir -p logs
./env/bin/python pipeline/p2_metrics.py --nii-dir data/nlst_serial/nifti --seg-dir data/seg/nlst_serial \
  --out results/p2/p2_paraspinal.csv --workers 60
