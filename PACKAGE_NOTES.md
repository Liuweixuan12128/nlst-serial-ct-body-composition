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


## Update on 2026-10-25

Source: the Code directory of final_JCSM_2026-10-24, supplied with the manuscript now being formatted for Frontiers in Endocrinology. This is a preparation snapshot, not evidence of journal submission or acceptance.

Compared with the previous repository revision, the only scientific source-code difference in the supplied package is figures/p2_figs_extra_rows.py. Its internal P2_Fig6 panel f (display Figure 7F) uses the observed upper confidence limit in the final-year window and power_grid.csv to report the largest compatible extra annual change for intramuscular fat and muscle volume. The grid uses onset one year before diagnosis and one extra normal annual rate as its calibration. The panel exports the upper z limit, the injected shift, the normal annual rate, the multiple and the corresponding extra annual change as source data. It no longer uses max_upper_bound_true_sd or the previous three-measure display. This change was supplied by the authors; the analysis was not rerun here.

README.md now maps the current 10 main figures and 8 supplementary figures to their internal output names. The six already-portable shell/SLURM scripts were retained instead of reintroducing the personal cluster paths from the supplied source package. All other analysis scripts and all historical plans/hash records match the existing repository byte for byte.

The journal-formatted Word manuscript, supplementary material, cover letter, STROBE checklist and independent figure files remain in the local submission folder. They are outside this code-only repository, as in the previous release. No credentials, participant data, review decisions or generated statistical outputs were added. No license was invented or changed.

Validation: Python files parsed without executing the analysis; changed plotting file compared byte for byte with the supplied revised source; all three historical plan hashes checked against their records; package file hashes regenerated. Numerical results, GPU segmentation, R models and figure rendering were not reproduced because the required data and results are absent from this snapshot.

## Supplementary-methods documentation update on 2026-10-25

The local submission supplement was revised as `paper2_Supplementary_Material_261026_revised.docx`. Its opening sections now describe the implemented study methods instead of reproducing project instructions, shell commands and internal review narratives. Supplementary Methods 1-4 retain their topic numbering. Analysis timing is summarised separately, and the descriptive thoracic slab-length findings are placed in supplementary results. The eight supplementary figures and 43 rows of result values were retained without recalculation.

`docs/paper2/SUPPLEMENTARY_METHODS_GUIDE.md` was added to connect the revised sections to the source code and historical plans. It records the distinction between the longitudinal and complete-component score samples, the separate normative and reliability models, annualised interval changes, and the approximate post hoc injection-recovery procedure. The counts and result values cited in the guide are from the supplied manuscript materials; they were not reproduced during this update. README.md links to this guide and states the outstanding interpretation check.

The body-size code defines negative log slab length, standardises it, and then reverses its sign again for the conditional logistic model, but labels the contrast as "per 1 SD shorter". This is documented as an unresolved discrepancy between the implemented predictor and its label. The supplement marks the issue for author review. The source model and reported numerical estimates remain unchanged. The supplement's two editorial replacements of "as published" with "main model" do not alter the labels emitted by the preserved statistical script; the guide makes this difference explicit.

All historical planning documents and their original hash records remain byte-identical to the preceding repository version. No statistical source file was changed for this supplementary-methods revision. The earlier supplied change to `figures/p2_figs_extra_rows.py` remains part of the same upload. The submission DOCX files and editorial notes remain outside the code repository.

Validation for this update covers the package checksums, all three historical plan hashes, Python syntax, local documentation links, the exact set of changes against the current remote revision, and archive contents. These checks do not establish numerical reproduction or resolve the slab-length contrast.
