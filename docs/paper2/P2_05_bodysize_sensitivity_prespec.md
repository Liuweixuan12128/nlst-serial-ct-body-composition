# Paper 2: body-size sensitivity analysis, declared before computation

Written 2026-10-18, after every result of `P2_00`, `P2_03` and `P2_04` was known. The SHA-256 of this file is recorded
in `results/p2/bodysize_prespec_sha256.txt` before `stats/p2_bodysize_sens.R` runs, and the script refuses to run if
the file has changed. Every result is reported. The analysis is not confirmatory.

## Why
Muscle volume depends on body size, and the public NLST tables contain no height or weight. The reference model for
the muscle-volume z-score adjusts for age, sex, smoking status and vendor only, so a lower muscle volume in cases could
reflect a smaller frame rather than less muscle. Reviewers of a muscle journal will ask for body-size normalisation.

## Proxy, and what was seen before this plan
The T6–T11 slab length (number of slices × slice spacing, `psm_nslices × slice_mm`, from the T11 to the T6 vertebral
centroid) measures the length of five thoracic segments, which scales with stature. Before this plan was written, the
following was looked at on each participant's first analysed scan (`results/p2/partA_long.csv`):
- men/women ratio of mean slab length 1.09; within-person coefficient of variation between rounds 0.9%, between-person
  8.1%;
- within-sex correlation of log slab length with log muscle volume 0.48, with log lung volume 0.24–0.37, with T10
  trabecular attenuation −0.01;
- mean slab length in cases versus controls: men 128.8 versus 128.5 mm, women 117.7 versus 117.7 mm.
No case–control model including slab length had been fitted.

## Analyses
Definitions as in `P2_03` (noise adjustment, reference model on controls' first scans, score S, OR per 1 SD lower,
conditional logistic regression on matched sets), original control sets and QC.
1. **Slab length itself.** OR per 1 SD shorter slab length (z against the same reference model as the other components:
   age, sex, smoking status, vendor), to show whether cases had a smaller frame.
2. **Muscle volume adjusted for body size.** The muscle-volume z is recomputed with log slab length added to its
   reference model. Refitted with this z:
   - muscle volume alone;
   - the primary score S;
   - the bone-plus-muscle score, alone, adjusted for the 15th-percentile lung density z, and adjusted for the
     %LAA−950 z (as defined in `P2_04` analysis 6).

## Reporting rule, set in advance
- Reported in one sentence in the Results and in full in Supplementary Table 1, whatever they show.
- If the muscle-volume OR adjusted for slab length has a CI that includes 1, the manuscript will say that the lower
  muscle volume in cases may reflect body size.
- The rule of `P2_04` (no bone or muscle signal independent of emphysema) stands regardless of this analysis.
