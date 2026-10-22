# Paper 2: sensitivity analyses prompted by internal review, declared before computation

Written 2026-10-15. Every result of `P2_00` (primary hypothesis) and `P2_03` (reserve score) was already known when
this plan was written. The SHA-256 of this file is recorded in `results/p2/rev_prespec_sha256.txt` before
`stats/p2_review_sens.R` runs, and the script refuses to run if the file has changed. Every result is reported. None of
these analyses is confirmatory. They test whether the reported estimates survive specific criticisms.

## Why
An internal review (2026-10-15) raised four points that the existing analyses did not answer.
1. **Control eligibility.** Controls had to be imaged in all three rounds, while a case needed only ≥2 rounds before
   diagnosis. For a case diagnosed before T2, a control therefore had to attend screening after the case's diagnosis.
2. **Visual quality control.** The visual review that excluded 14 scans was performed by an AI assistant that could
   see case/control labels. It has not yet been confirmed by a radiologist.
3. **Discrimination of lung density.** The concordance of lung density alone was not reported, so it is unknown
   whether the composite adds discrimination beyond it.
4. **Lung density measurement.** The 15th-percentile density was not adjusted for inspiration level (lung volume).
   The standard emphysema index (%LAA below −950 HU) was not reported.
The review also asked how strong unmeasured confounding would have to be (E-value).

## Analyses
All use the definitions of `P2_03` (noise adjustment, control reference model on baseline scans, score S, OR per 1 SD
lower, conditional logistic regression on matched sets) unless stated.

1. **Risk-set controls.**
   - Controls were re-drawn by `pipeline/p2_resample_riskset.py` (seed 20261015). This was done before any outcome
     analysis of the new controls.
   - Eligibility: never diagnosed with lung cancer, still free of it on the case's diagnosis day, and imaged at the
     rounds the case had before diagnosis. Matching factors are unchanged.
   - Within each case's risk set, a hypergeometric number of the original controls was kept, and the rest were drawn
     from the newly eligible participants. Result: 45 of 2,175 controls replaced, all in sets whose case had two
     pre-diagnostic rounds.
   - New scans go through the same segmentation, measurement and automatic QC (slab-length and prone rules). They
     have had no visual review.
   - `stats/p2_partA.R` is run with these control sets (`P2_CTRL=_riskset`, `P2_BOOT=0`) into
     `results/p2/sens_riskset`.
   - Reported from that run: the primary-hypothesis estimates (interval z of intramuscular fat by window) and, from the
     score, the primary OR, the components, the musculoskeletal score adjusted for lung density, and concordance.
2. **No visual-QC exclusions.** `stats/p2_partA.R` with the visual-review exclusions removed (automatic exclusions
   kept), into `results/p2/sens_novisualqc`. Reported: the same estimates as analysis 1.
3. **Sets whose case had three pre-diagnostic rounds.** Here the control eligibility of the original selection was
   already correct. Primary score model in the original data, restricted to these sets.
4. **Discrimination.** Matched-set concordance of
   - S,
   - lung density z alone,
   - the musculoskeletal score alone.
   For each, a matched-set bootstrap 95% CI (1,000 replicates, seed 20261015), plus the paired bootstrap difference
   between S and lung density alone.
5. **Lung density adjusted for lung volume.** z_L is recomputed with log lung volume added to the control reference
   model. The primary score model is refitted with this z_L, and the lung-only model with it.
6. **%LAA below −950 HU** (`lung_laa950_pct`, log(x + 0.1) transformed, oriented so that higher = more reserve).
   It replaces Perc15 in z_L, and the lung-only and primary score models are refitted.
7. **E-values** (VanderWeele and Ding) for the primary OR and for the musculoskeletal OR adjusted for lung density, for
   the point estimate and the confidence limit nearer 1. The OR is treated as a risk ratio, because cumulative
   incidence in the CT arm is about 4%.

## Reporting rule, set in advance
- The original analysis remains the main analysis. Analyses 1–7 are reported in the Results (one sentence each) and
  in full in a supplementary table, whatever they show.
- If the primary OR in analysis 1 or 2 has a CI that includes 1, the abstract will say so.
- If the musculoskeletal OR adjusted for lung density has a CI that includes 1 in analysis 1, 2, 5 or 6, the
  manuscript will not describe a musculoskeletal signal independent of emphysema.
