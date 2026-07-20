# autoc_hte

Heterogeneous treatment effects (HTE) of intensive BP control on eGFR slope,
evaluated with the **AUTOC** (Area Under the TOC curve) framework across SPRINT,
ACCORD, AASK, and MDRD.

**Primary estimand: chronic eGFR slope** — the post–acute-dip slope (`time2` in a
piecewise-linear spline with the acute/chronic knot fixed at 4 months). Total
slope (the whole-trajectory slope) was the prior main analysis and is reported
here only as a secondary comparator.

Patients are ranked by a baseline covariate (baseline eGFR, **UACR**, or KFRE);
at each percentile cutoff the LMM is refit on that subgroup to read off the
subgroup CATE, and the **AUTOC score** is the area under the (subgroup CATE −
overall CATE) curve — larger |AUTOC| ⇒ stronger effect modification by that
covariate. **UACR was found to be the strongest HTE marker** (step 3), and is
the covariate carried into the threshold work (steps 4–5b).

## Pipeline & run order

| Step | Script | Purpose | Output |
|---|---|---|---|
| — | `R/00_functions.R` | Shared `calculate_slopes()` (chronic-slope LMM) + `auc_trapezoid()` | sourced by others |
| 1 | `R/01_build_analysis_data.R` | Pool the trials + baseline covariates into one long dataset | `data/derived/analysis_long.xlsx` |
| 2 | `R/02_autoc_chronic.R` | Single-pass AUTOC point estimates (eGFR / UACR / KFRE), pooled + per trial | `results/autoc_chronic_point.csv`, `results/chronic_ate_point.csv` |
| 3 | `R/03_autoc_bootstrap.R` | Patient-bootstrap AUTOC CIs, pairwise covariate comparisons, winner tally, UACR thresholds & subgroup CATEs (heavy → Sherlock array) | `frag_*` fragments |
| 3b | `R/03b_combine_bootstrap.R` | Stitch step-3 fragments → summaries + 95% CI curve ribbon | `results/autoc_bootstrap_results_<COHORT>.xlsx`, `results/autoc_curves_for_plot_<COHORT>.csv` |
| 4 | `R/04_cv_point.R` | **Threshold Part 1**: full-data pooled UACR thresholds (chronic CATE>1.0, total CATE>0.75), no bootstrap | `results/cv_point.csv`, `results/cv_point_curve.csv` |
| 5 | `R/05_autoc_cv.R` | **Threshold Part 2**: 1000× Monte-Carlo CV — 60/40 trial-stratified split, find UACR threshold in TRAIN, estimate held-out above/below effect in TEST (heavy → Sherlock array) | `frag_cv_*` fragments |
| 5b | `R/05b_combine_cv.R` | **Threshold Part 3**: stitch CV fragments → held-out effect 95% CIs, threshold stability, success rate, and TRAIN−TEST optimism/bias | `autoc_cv_results.xlsx`, `autoc_cv_raw.csv` |

`sherlock/` holds the SLURM `.sbatch` wrappers for the heavy array steps
(`run_autoc_bootstrap_array.sbatch` for step 3, `run_autoc_cv_array.sbatch` for
step 5) and the submit/combine recipes.

## Analyses done so far

1. **Chronic-slope LMM & overall ATE** (steps 0–1) — piecewise spline, 4-month
   knot, per-patient random slopes; pooled and per-trial chronic-slope ATE of
   intensive vs standard BP control.
2. **AUTOC HTE screen** (steps 2–3b) — which baseline covariate (eGFR, UACR,
   KFRE) best captures effect modification on the chronic (and total) slope,
   with patient-bootstrap CIs and pairwise comparisons. **Result: UACR wins.**
3. **UACR threshold discovery & cross-validation** (steps 4–5b) —
   - *Part 1 (04):* on the full pooled cohort, the UACR value where the chronic
     CATE crosses +1.0 and where the total CATE crosses +0.75 (single point
     estimate + in-sample above/below effects).
   - *Part 2 (05):* Monte-Carlo cross-validation — 1000 independent 60/40
     patient-level splits, stratified by trial, rediscovering the threshold in
     each 60% training fold and testing it on the disjoint 40%.
   - *Part 3 (05b):* honest held-out treatment effect above vs below the
     threshold (95% CI) for both slopes, plus the reliability diagnostics
     (crossing/success rate, threshold distribution) and the **optimism/bias**
     of the AUTOC process (in-sample TRAIN effect minus held-out TEST effect).

### Cross-validation design notes (steps 5 / 5b)

- **Split is by patient `id`, stratified by trial, without replacement** — a
  patient's whole trajectory stays on one side (no leakage), and every split
  keeps the pooled SPRINT/ACCORD/AASK/MDRD mix.
- **Monte-Carlo CV** (repeated random 60/40 holdout) was chosen over .632
  bootstrap-validation because the disjoint, no-replacement test fold most
  directly answers "does the discovered threshold hold up on genuinely
  different patients?" without duplicated patients contaminating the LMM's
  random-effects structure.
- **Failed iterations are kept in the denominator.** An iteration where the
  target CATE is never crossed in TRAIN (or where a TEST subgroup fit fails to
  converge) is excluded from the effect summaries but still counted toward the
  success rate — the success rate is itself a primary result, not discarded.
- **Bias check.** The paired TRAIN−TEST effect gap per iteration estimates the
  optimism of the discovery process; the tiers full-data point (04) →
  in-sample TRAIN → held-out TEST show the shrinkage from optimistic to honest.

## Data

All inputs live in `data/raw/` and are **git-ignored** (IRB — patient-derived,
never pushed). See [data/README.md](data/README.md) for provenance. Scripts read
via the relative `data/` path; no `setwd()` to absolute machine paths. On
Sherlock the heavy steps read `analysis_long.xlsx` from `$HOME` (`AUTOC_DATA`).

## Provenance

Migrated 2026-06-25 from the original working folder `Project 2_ HTE Analysis/`.
Current code came from `New LMM code_AUTOC/` and `Sherlock scripts and results/`;
superseded datasets and side analyses were left behind. Threshold
cross-validation (steps 4–5b) added 2026-07-20.
