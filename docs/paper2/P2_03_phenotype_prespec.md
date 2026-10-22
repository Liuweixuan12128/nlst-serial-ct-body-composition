# Paper 2 — CT host-phenotype score: analysis plan declared before computation

Written 2026-10-11, before any of the analyses below were run. The SHA-256 of this file is recorded in
`results/p2/pheno_prespec_sha256.txt` before `stats/p2_phenotype.R` is executed; the script refuses to run if the
file has changed since. The analysis is run once and every result is reported.

## Status and disclosure (to be stated in the manuscript)
- The pre-specified primary hypothesis (accelerated muscle ageing before diagnosis, `P2_00` section 4) was null.
- The pre-specified trajectory model (`group × tau`) contained a group main effect. It showed lower bone density,
  lower muscle volume and lower lung density in cases at every time point (`P2_00` section 8.4).
- This plan was written **after** those results were known. The composite therefore summarises component differences
  that had already been observed; it is not an independent confirmation of them. What it adds, and what was not yet
  known, is: the size of the joint association, its gradient across quartiles, whether the musculoskeletal part is
  independent of emphysema, whether it replicates in the two independently sampled case series, whether it is present
  ≥4 years before diagnosis, and whether the choice of components matters.

## Population
Exactly the analysis population of `stats/p2_partA.R` (`results/p2/partA_long.csv`): QC-passing scans; for cases, scans
on or before the diagnosis day; participants with ≥2 QC-passing scans; 714 cases and 2,150 matched controls in
725 matched sets. Each participant's **baseline scan** = first QC-passing scan (smallest `t`).

## Components (oriented so that higher = more reserve)
| code | measure | variable |
|---|---|---|
| B | T10 trabecular bone radiodensity (HU) | `T10_trab_hu` |
| M | paraspinal muscle volume, T6–T11 (100 × log cm³) | `psm_lvol` |
| L | lung density, 15th percentile (HU; higher = less emphysema) | `lung_perc15_hu` |

For each component:
1. Noise adjustment as in the main analysis: y − γ·`noise_c`, with γ = the normative model's noise coefficient
   (`results/p2/partA_normative_model.csv`, term `noise_c`).
2. Reference model fitted on **controls' baseline scans only**: `lm(y ~ age0c + sex + smoke + vendor0)`.
3. z = (y − prediction) / residual SD of that control model, for cases and controls alike.

**Score S** = mean of z_B, z_M, z_L (all three required). It is then divided by its SD in controls, so that
1 unit = 1 control SD. **Primary exposure: S reversed (−S), i.e. per 1 SD lower reserve.**

## Primary analysis
Conditional logistic regression on matched sets: `clogit(case ~ I(-S) + strata(set))`; OR per 1-SD lower score,
95% CI, Wald P. Two-sided α = 0.05. This is the only test of this plan treated as confirmatory.

## Supporting analyses (estimates with 95% CI; no multiplicity adjustment; labelled as supporting)
1. **Quartiles**: cut-points = quartiles of S among controls. Q4 (highest reserve) = reference. ORs for Q3, Q2, Q1;
   trend test with quartile index (1–4) as a continuous term.
2. **Individual components** in the same clogit framework (OR per 1-SD lower z), for comparison with the composite.
3. **Musculoskeletal score** S_MSK = mean(z_B, z_M), standardised in controls: (a) alone; (b) adjusted for z_L. This tests
   whether bone and muscle carry information beyond emphysema.
4. **Component-choice sensitivity**: a five-measure score using all five pre-specified measures of the main analysis
   (z_B, z_M, z_L, intramuscular fat % reversed, lean radiodensity `psm_hu`), built the same way.
5. **Internal replication**: primary model fitted separately in the resected series (`p2`) and the unresected series
   (`p2b`). These are independently sampled case series, and no control appears in both. Heterogeneity is tested by a
   likelihood-ratio test of the S × series interaction.
6. **Lead time**: matched sets split by the years from the case's baseline scan to diagnosis: <2, 2–<4, ≥4 years.
   The primary model is fitted in each stratum, with a likelihood-ratio test of the S × stratum interaction.
7. **Latest pre-diagnostic scan**: the primary model repeated using each participant's last QC-passing scan (cases:
   last scan on or before diagnosis), with the control reference model refitted on controls' last scans.
8. **Exploratory subgroups** (S × subgroup interaction, LR test): histology (adenocarcinoma, squamous, small cell,
   other), stage (I–II vs III–IV), detection (screen-detected vs interval/post-screening), sex, smoking status.
9. **Discrimination**: the matched-set concordance of S (share of case–control pairs within a set in which the case has
   the lower score; ties count ½), with a matched-set bootstrap 95% CI (1,000 replicates, seed 20261011).

## Outputs
`results/p2/pheno_scores.csv` (per participant), `pheno_estimates.csv` (every model above), `pheno_quartiles.csv`,
`pheno_summary.txt`.

## Interpretation, set in advance
- OR > 1 with a CI above 1 → a lower combined bone, muscle and lung reserve on screening CT marks people who are
  later diagnosed with lung cancer. This is a marker of risk, not evidence that low reserve causes cancer.
- Residual confounding by cumulative smoking exposure (pack-years are not in the public tables; matching was on
  current versus former smoking only) is the main alternative explanation. It is stated as such and is not tested
  here.
- If the musculoskeletal score is null after adjustment for z_L, the result is reported as an emphysema signal and
  the manuscript does not claim a musculoskeletal phenotype.
