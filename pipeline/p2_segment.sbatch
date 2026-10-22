#!/bin/bash
#SBATCH --job-name=p2seg
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --array=0-5
#SBATCH --time=12:00:00
#SBATCH --output=logs/p2_seg.%A_%a.out
# Paper 2 segmentation: 6 GPUs x 3 workers, all sharing one folder through lock files (p2_segment_worker.py).
P2_PROJECT_ROOT="${P2_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$P2_PROJECT_ROOT"
mkdir -p logs
for w in 0 1 2; do
  ./env/bin/python pipeline/p2_segment_worker.py --nii-dir data/nlst_serial/nifti --seg-dir data/seg/nlst_serial \
    --done-flag ${DONE_FLAG:-logs/p2_fetch.log} > logs/p2_seg${LOG_TAG}.${SLURM_ARRAY_TASK_ID}.${w}.log 2>&1 &
  sleep 20
done
wait
