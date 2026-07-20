# Data — provenance

**IRB / privacy:** every file in `raw/` is patient-derived and is git-ignored.
Nothing in this folder is pushed to GitHub. Keep these local only.

These are the **only** datasets the current chronic-slope AUTOC analysis uses,
pulled from the original `HTE Datasets/` working folder. Copied 2026-06-25.

## raw/ — slope & AUTOC inputs (long format, one row per eGFR visit)

| File | Source (original location) | Cohort | Key columns | n rows |
|---|---|---|---|---|
| `combined_long.xlsx` | `Excel Datasets/combined_long.xlsx` (2025-06-18) | SPRINT + ACCORD | `id, arm, egfr, days, uacr, kfrs, egfr0` | ~20,789 |
| `aask_dataset.xlsx` | `Excel Datasets/aask_dataset.xlsx` (2025-09-05) | AASK | `id, arm, egfr, days, uacr, kfrs, egfr0` | ~7,303 |
| `mdrd_dataset.xlsx` | `Excel Datasets/mdrd_dataset.xlsx` (2025-09-16) | MDRD | `id, arm, egfr, days, uacr, kfrs, egfr0` | ~8,114 |

Notes:
- In `combined_long.xlsx`, trial is the first character of `id` (`S` = SPRINT, `A` = ACCORD).
- `aask_dataset.xlsx` / `mdrd_dataset.xlsx` are the *self-contained* successors to the
  older `aask_0630.xlsx` + `aask_long(4).xlsx` (and MDRD equivalents): UACR / eGFR0 / KFRS
  are already merged in, so no separate baseline-UACR join is needed. The `_0630` and
  `_long(4)` files are **superseded** and intentionally not carried over.

## raw/ — outcomes (one row per participant; for HRs & CIF figures)

| File | Source | Cohort |
|---|---|---|
| `sprint_events.xlsx` | `HTE Datasets/sprint_events.xlsx` | SPRINT |
| `accord_events.xlsx` | `HTE Datasets/accord_events.xlsx` | ACCORD |
| `aask_events.xlsx`   | `HTE Datasets/aask_events.xlsx`   | AASK |
| `mdrd_events.xlsx`   | `HTE Datasets/mdrd_events.xlsx`   | MDRD |

Each holds time-to-event pairs (`fu_*`, event flags) for death, CVD, KRT, etc.
Note `mdrd_events.xlsx` uses `ID` (uppercase) — lowercase it on read.

## raw/ — Table 1 covariates

| File | Source | Contents |
|---|---|---|
| `all_jul2.xlsx` | `Excel Datasets/all_jul2.xlsx` | Baseline covariates (age, sex, BMI, comorbidities, labs) for all 4 cohorts |
