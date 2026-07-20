# =============================================================================
# 00_functions.R  --  Shared functions for the chronic-slope AUTOC analysis.
#
# This file defines functions only; it does not run an analysis by itself.
# Every analysis script loads it with:  source("R/00_functions.R")
#
# Two functions are provided:
#   calculate_slopes() : fit the piecewise LMM, return the CHRONIC-slope ATE
#   auc_trapezoid()    : trapezoidal area-under-the-curve (the AUTOC score)
# =============================================================================

suppressPackageStartupMessages({
  library(nlme)
  library(lspline)
  library(dplyr)
  library(tibble)
})

# Acute/chronic boundary (the spline knot), in MONTHS, used uniformly for all
# cohorts including SPRINT. lspline() takes the knot in days, so we use knot*31.
KNOT_MONTHS <- 4

# -----------------------------------------------------------------------------
# calculate_slopes()
#   Fit:   egfr ~ time1*arm + time2*arm,  random = ~1 + time1 + time2 | id
#   where time1 = pre-knot ("acute") segment, time2 = post-knot ("chronic").
#
#   Returns per-arm CHRONIC (time2) slopes per year and the chronic-slope ATE
#   (treatment - control). Total slope is intentionally not reported: this
#   pipeline is chronic-slope only.
# -----------------------------------------------------------------------------
calculate_slopes <- function(dataset, knot = KNOT_MONTHS) {

  dataset$time1 <- lspline(dataset$days, knot * 31)[, 1]
  dataset$time2 <- lspline(dataset$days, knot * 31)[, 2]

  fit <- lme(egfr ~ time1 * arm + time2 * arm,
             data    = dataset,
             random  = ~ 1 + time1 + time2 | id,
             control = lmeControl(opt = "optim"))

  fixed  <- fixef(fit)
  random <- ranef(fit)

  # Find the time2:arm interaction term robustly: R may label it "time2:arm"
  # or "arm:time2" depending on term ordering, so match on both pieces.
  nm        <- names(fixed)
  t2arm_nm  <- nm[grepl("time2", nm) & grepl("arm", nm)]
  beta_t2   <- fixed["time2"]
  beta_t2a  <- fixed[t2arm_nm]

  # Patient-specific chronic (time2) slopes, expressed per year (days * 365).
  patient_slopes <- tibble(id = rownames(random)) %>%
    mutate(
      chronic_noarm = (beta_t2 + random[id, "time2"]) * 365,
      chronic_arm   = (beta_t2 + beta_t2a + random[id, "time2"]) * 365,
      chronic_diff  = chronic_arm - chronic_noarm
    )

  population_effects <- patient_slopes %>%
    summarize(
      population_chronicslope_noarm = mean(chronic_noarm),
      population_chronicslope_arm   = mean(chronic_arm),
      population_chronicslope_diff  = mean(chronic_diff)
    )

  list(
    model_summary      = summary(fit),
    patient_slopes     = patient_slopes,
    population_effects = population_effects
  )
}

# -----------------------------------------------------------------------------
# auc_trapezoid()
#   Area under the curve (x, y) by the trapezoidal rule. Used to integrate the
#   (subgroup ATE - overall ATE) curve across percentiles into the AUTOC score.
# -----------------------------------------------------------------------------
auc_trapezoid <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
