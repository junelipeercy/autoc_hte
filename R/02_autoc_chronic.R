# =============================================================================
# 02_autoc_chronic.R  --  AUTOC on the CHRONIC slope, single run (no bootstrap).
#
# Question: is there heterogeneity of treatment effect (HTE) on the chronic
# eGFR slope, ranked by each of three baseline covariates (eGFR, UACR, KFRE)?
#
# For each cohort (pooled, then each trial separately) and each covariate:
#   - rank patients by the baseline covariate
#   - at each percentile cutoff, refit the LMM on that subgroup and read off
#     the chronic-slope ATE (treatment - control)
#   - AUTOC score = area under the (subgroup ATE - overall ATE) curve.
#     |AUTOC| near 0  -> little HTE;  larger -> stronger HTE on that covariate.
#
# No subgroup split, no bootstrap, no CIs here -- one pass to see the signal.
#
# Run 01_build_analysis_data.R first to create data/derived/analysis_long.xlsx.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})
source("R/00_functions.R")

KNOT  <- KNOT_MONTHS   # 4 months, uniform for all cohorts
N_PCT <- 90            # number of percentile cutoffs to integrate over

# -----------------------------------------------------------------------------
# Overall chronic-slope ATE for a cohort (whole-cohort fit, no subgrouping).
# -----------------------------------------------------------------------------
overall_chronic_ate <- function(df, knot = KNOT) {
  eff <- calculate_slopes(df, knot)$population_effects
  c(chronic_control   = eff$population_chronicslope_noarm,
    chronic_treatment = eff$population_chronicslope_arm,
    chronic_ate       = eff$population_chronicslope_diff)
}

# -----------------------------------------------------------------------------
# AUTOC score for one cohort + one covariate.
#   biomarker: "egfr" (rank low->high risk = low eGFR), "uacr" or "kfre"
#              (rank high value = high risk).
#   Returns the AUTOC score and the underlying curve for plotting later.
# -----------------------------------------------------------------------------
compute_autoc <- function(df, biomarker, knot = KNOT, n_pct = N_PCT) {

  # baseline covariate value is constant per patient -> one row per id
  base_vals <- df %>%
    distinct(id, .keep_all = TRUE) %>%
    transmute(val = switch(biomarker,
                           egfr = egfr0,
                           uacr = uacr,
                           kfre = kfrs))
  vals <- sort(unique(base_vals$val[!is.na(base_vals$val)]))
  if (length(vals) < n_pct) return(list(auc = NA_real_, curve = NULL))

  pct <- quantile(vals, probs = seq(0.01, 1, by = 0.01), na.rm = TRUE)
  # eGFR: lower = higher risk, so walk percentiles top-down; UACR/KFRE bottom-up
  cutoffs <- if (biomarker == "egfr") sort(pct, decreasing = TRUE)[1:n_pct]
             else                      sort(pct)[1:n_pct]

  keep <- NULL
  for (i in cutoffs) {
    sub <- switch(biomarker,
                  egfr = df[df$egfr0 <= i, ],
                  uacr = df[df$uacr  >= i, ],
                  kfre = df[df$kfrs  >= i, ])
    sub <- sub %>% filter(!is.na(egfr), !is.na(days), !is.na(arm), !is.na(id))

    ate <- tryCatch(
      calculate_slopes(sub, knot)$population_effects$population_chronicslope_diff,
      error = function(e) NA_real_)
    keep <- rbind(keep, c(cutoff = i, ate = ate))
  }

  # AUTOC = area under (subgroup ATE - overall ATE). keep[1,] is the broadest
  # subgroup (~ whole cohort), so keep[1,"ate"] approximates the overall ATE.
  auc <- auc_trapezoid(n_pct:1, keep[, "ate"] - keep[1, "ate"])
  list(auc = auc, curve = as.data.frame(keep))
}

# -----------------------------------------------------------------------------
# Run everything: pooled cohort + each trial.
# -----------------------------------------------------------------------------
analysis_long <- read_excel("data/derived/analysis_long.xlsx")

cohorts <- list(
  Pooled = analysis_long,
  SPRINT = filter(analysis_long, trial == "SPRINT"),
  ACCORD = filter(analysis_long, trial == "ACCORD"),
  AASK   = filter(analysis_long, trial == "AASK"),
  MDRD   = filter(analysis_long, trial == "MDRD")
)

# Percentile depth per cohort. The deepest subgroups still converge at 90, but
# the smaller trials (ACCORD, MDRD) are capped at 85 so the extreme-percentile
# ATEs are less noisy. NOTE: AUTOC magnitude scales with this depth, so scores
# are directly comparable ACROSS COVARIATES within a cohort (same depth), but
# not across cohorts that use different depths.
cohort_npct <- c(Pooled = 90, SPRINT = 90, ACCORD = 85, AASK = 90, MDRD = 85)

autoc_rows  <- list()
ate_rows    <- list()
curves      <- list()   # kept in memory for later plotting

for (cname in names(cohorts)) {
  df    <- cohorts[[cname]]
  n_pct <- cohort_npct[[cname]]
  cat("\n=====================  ", cname,
      "  (n =", n_distinct(df$id), "patients, depth =", n_pct,
      "pct)  =====================\n")

  ate <- overall_chronic_ate(df)
  cat(sprintf("Overall chronic ATE: %.3f  (control %.3f, treatment %.3f) mL/min/1.73m2/yr\n",
              ate["chronic_ate"], ate["chronic_control"], ate["chronic_treatment"]))
  ate_rows[[cname]] <- data.frame(cohort = cname, t(ate))

  for (bm in c("egfr", "uacr", "kfre")) {
    res <- compute_autoc(df, bm, n_pct = n_pct)
    cat(sprintf("  AUTOC[%-4s] = %s\n", bm,
                ifelse(is.na(res$auc), "NA (too few patients)",
                       sprintf("%.3f", res$auc))))
    autoc_rows[[paste(cname, bm)]] <-
      data.frame(cohort = cname, covariate = bm, n_pct = n_pct, autoc = res$auc)
    curves[[paste(cname, bm)]] <- res$curve
  }
}

autoc_table <- do.call(rbind, autoc_rows); rownames(autoc_table) <- NULL
ate_table   <- do.call(rbind, ate_rows);   rownames(ate_table)   <- NULL

cat("\n\n================  AUTOC scores (chronic slope)  ================\n")
print(autoc_table, row.names = FALSE)
cat("\n================  Overall chronic-slope ATEs  ================\n")
print(ate_table, row.names = FALSE)

# Save point estimates (curves stay in memory as `curves` for plotting later)
dir.create("results", showWarnings = FALSE)
write.csv(autoc_table, "results/autoc_chronic_point.csv", row.names = FALSE)
write.csv(ate_table,   "results/chronic_ate_point.csv",   row.names = FALSE)
cat("\nWrote results/autoc_chronic_point.csv and results/chronic_ate_point.csv\n")
