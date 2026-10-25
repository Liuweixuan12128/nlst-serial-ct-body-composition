# Guide to the revised supplementary methods

Prepared on 2026-10-25 for the local submission supplement `paper2_Supplementary_Material_261026_revised.docx`. The filename is a local revision identifier. This guide maps the supplement to the code and records clarifications; it does not replace the original analysis plans or establish a new prospective plan. The submission documents are maintained outside this code package.

## Section-to-source mapping

| Revised supplement section | Implementation | Historical plan |
|---|---|---|
| Methods 1: CT reserve score | [p2_phenotype.R](../../stats/p2_phenotype.R) | [P2_03_phenotype_prespec.md](P2_03_phenotype_prespec.md) |
| Methods 2: Additional sensitivity analyses | [p2_review_sens.R](../../stats/p2_review_sens.R), [p2_resample_riskset.py](../../pipeline/p2_resample_riskset.py), and alternative-input runs of [p2_partA.R](../../stats/p2_partA.R) | [P2_04_review_sensitivity_prespec.md](P2_04_review_sensitivity_prespec.md) |
| Methods 3: Image processing and longitudinal analyses | [p2_select_serial.py](../../pipeline/p2_select_serial.py), [p2_segment_worker.py](../../pipeline/p2_segment_worker.py), [p2_metrics.py](../../pipeline/p2_metrics.py), [p2_qc.py](../../pipeline/p2_qc.py), [p2_partA.R](../../stats/p2_partA.R), [p2_map_stats.py](../../stats/p2_map_stats.py), and [p2_power.py](../../stats/p2_power.py) | Original `P2_00` and `P2_01` study documents in this directory; the later exploratory additions are disclosed separately |
| Methods 4: Body-size sensitivity analysis | [p2_bodysize_sens.R](../../stats/p2_bodysize_sens.R) | [P2_05_bodysize_sensitivity_prespec.md](P2_05_bodysize_sensitivity_prespec.md) |
| Supplementary Table 1: Additional sensitivity analyses | Outputs of [p2_review_sens.R](../../stats/p2_review_sens.R) and [p2_bodysize_sens.R](../../stats/p2_bodysize_sens.R) | Methods 2 and 4 above |

The existing README gives the execution order and maps the 10 main figures and eight supplementary figures to their internal output names. Methods 3.8 describes the post hoc calibration used in main Figure 7F (internal `P2_Fig6`, panel f); its plotting implementation is [p2_figs_extra_rows.py](../../figures/p2_figs_extra_rows.py).

## Clarifications in the revised supplement

- **Analysis populations.** The supplied manuscript reports 714 cases and 2,150 controls in the longitudinal population. The primary score association model requires complete components and a matched set containing a case and at least one eligible control; it reports 714 cases and 2,115 controls in 714 matched sets. The historical score plan's reference to 725 sets is retained only in the unchanged historical record. These sample counts were not recomputed here.
- **Normative ageing and reliability.** In `p2_partA.R`, the normative fixed-effects model includes interactions of time with age, sex and smoking status. The separate reliability model omits those interactions. Normative coefficient intervals use the Wald method; slope-reliability intervals use participant bootstrap resampling. The revised prose distinguishes the two models.
- **Interval and trajectory models.** Interval changes are annualised before control-reference standardisation. The time-to-diagnosis model includes participant random intercepts and slopes and a matched-set random intercept. Its group term describes the case-control contrast at the diagnosis time origin.
- **Alternative controls.** The historical filename `p2_resample_riskset.py` is preserved. Its eligible pool remains restricted to participants never diagnosed with lung cancer; the analysis relaxes screening attendance requirements and is not unrestricted incidence-density sampling.
- **Injection-recovery analysis.** This analysis was added after the primary findings were known. It uses matched-set bootstrap calculations of interval contrasts and does not refit the primary mixed model in every replicate; detection probabilities are approximate. It must not be described as planned before any results were seen.
- **Editorial changes.** Operational commands, output paths and project instructions were removed from the submission supplement, while methodological definitions and analysis timing were retained. Thoracic slab-length descriptive findings were moved to supplementary results. The eight figures, their captions and 43 rows of numerical results were retained. Two table labels now read "main model" in place of "as published"; the preserved body-size script still emits the older labels.

## Slab-length direction check

The supplied body-size script contains these transformations:

```r
y_slab = -log(slab_mm)
z_slab = zfit("y_slab")
nz_slab = -z_slab
```

It then fits `case ~ nz_slab + strata(set)` but labels the estimate `T6-T11 slab length, per 1 SD shorter`. Because the standardised negative-log residual is negated again, the implemented predictor increases with residual log slab length. This conflicts with the stated shorter-length direction.

The source code and the actual outputs used to produce the reported estimate must be reconciled before choosing a correction. The revised supplement highlights the method and table label for that check and preserves the reported OR of 1.02 (95% CI 0.93-1.11; P = 0.73). This update does not invert that estimate, change the model or claim that the label has been resolved. Required per-scan and intermediate analysis files are absent from this source snapshot.

## Timing and unchanged historical records

The primary protocol was not registered. The reserve-score plan was dated 11 October 2026 after component differences had been examined on 9 October. The further sensitivity plans were dated 15 and 18 October after preceding results were known, and descriptive slab-length relationships had been examined before the body-size plan. The score and subsequent sensitivity analyses are therefore exploratory; the two case series provide consistency comparisons rather than independent confirmation of a newly defined score.

All original plan files and the three `results/p2/*_prespec_sha256.txt` records are preserved byte for byte. Their locally recorded dates and SHA-256 values identify historical files; they do not establish independent preregistration. No study model was rerun during this editorial and repository update.
