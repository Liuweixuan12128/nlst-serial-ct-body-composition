#!/bin/bash
#SBATCH --job-name=p2met
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --array=0-3
#SBATCH --time=12:00:00
#SBATCH --output=logs/p2_metrics.%A_%a.out
# Paper 2 measurement, sharded over nodes. PASS=n names the output set; SKIP=<csv> skips scans measured in earlier passes.
# Merge afterwards with pipeline/p2_metrics_merge.py.
P2_PROJECT_ROOT="${P2_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$P2_PROJECT_ROOT"
mkdir -p logs
mkdir -p results/p2/shards
./env/bin/python pipeline/p2_metrics.py --nii-dir data/nlst_serial/nifti --seg-dir data/seg/nlst_serial \
  --out results/p2/shards/pass${PASS}_${SLURM_ARRAY_TASK_ID}.csv --workers 44 \
  --shard ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} ${SKIP:+--skip $SKIP}
