# Package notes

Prepared: 2026-10-22 (Asia/Shanghai).
Source: the Code directory in the 2026-10-19 submission snapshot.

## Selection

Included: all files from pipeline/, stats/, figures/, docs/ and results/ in that Code directory. The original README was replaced with a repository-specific overview and the available reproduction prerequisites.
Excluded: manuscript/ (document builder, manuscript/supplement/cover-letter/checklist sources, bibliography helper and Word style template); the external submission documents and generated figure files; the original submission ZIP; review working files.

## Changes to existing source files

The following six files were changed only to remove the original personal cluster working-directory and log-output paths. SLURM examples now use P2_PROJECT_ROOT or SLURM_SUBMIT_DIR; the statistics wrapper resolves its root from its own location. The examples create logs/ after entering the root. For SLURM, create logs/ before submission as explained in README.md. Resource requests, script arguments and model code were retained.

- `pipeline/p2_fov.sbatch`
- `pipeline/p2_maps.sbatch`
- `pipeline/p2_metrics.sbatch`
- `pipeline/p2_metrics_shard.sbatch`
- `pipeline/p2_segment.sbatch`
- `stats/p2_run_all.sh`

All other included source files and all historical plan/hash files remain byte-identical to the supplied code snapshot. The .gitattributes policy prevents automatic newline conversion of the hashed plans.

## Inputs absent from the supplied snapshot

- Raw CT images, NIfTI conversions and segmentation masks.
- The public clinical-data cache and exported clinical/selection tables.
- Per-scan measurements, including results/p2/p2_paraspinal.csv.
- The original visual-review decisions in results/p2/p2_qc_review.csv and a completed blinded radiologist review.
- Generated Part A, score, sensitivity, power and source-data outputs.
- Image/render assets used in the final figures.

These absent files were not replaced with dummy or simulated data. This is a source-code package, not a complete numerical-results archive. The original downloader's reuse cache and the final QC decisions require preparation before an end-to-end reconstruction.

## Verification performed

- Python source parsed for syntax without executing the analysis.
- Unchanged copied files compared with the source using SHA-256.
- All three analysis-plan SHA-256 values matched their original records.
- Package checked for excluded submission files and original personal cluster paths.
- ZIP contents checked against the prepared folder.

GPU processing, R models and end-to-end numerical reproduction were not executed as part of file packaging.
