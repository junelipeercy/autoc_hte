# =============================================================================
# 04_cv_point.R  --  PART 1 of the threshold-validation analysis.
#
# ONE pass on the FULL pooled dataset (no bootstrap, no split). It:
#   - builds the UACR-ranked CATE curve (chronic + total slope) on everyone,
#   - reports the pooled AUTOC score for UACR (both slopes),
#   - reads off the two clinically-motivated UACR thresholds:
#         chronic slope : the UACR value where chronic CATE first crosses +1.00
#         total   slope : the UACR value where total   CATE first crosses +0.75,
#   - splits the cohort at each threshold and reports the above/below subgroup
#     treatment effects on the WHOLE sample (the "optimistic" in-sample number
#     that the cross-validation in 05_autoc_cv.R is meant to honestly correct).
#
# This is light -- run it locally / on an sh_dev node, not the array.
# Run 01_build_analysis_data.R first to create data/derived/analysis_long.xlsx.
#
# Env vars (optional):
#   AUTOC_DATA   path to analysis_long.xlsx (default data/derived/analysis_long.xlsx)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(nlme); library(lspline); library(dplyr)
})

KNOT           <- 4
N_PCT          <- 90       # percentile depth on the pooled cohort (matches step 3)
TARGET_CHRONIC <- 1.00     # UACR threshold: chronic CATE crosses +1.00
TARGET_TOTAL   <- 0.75     # UACR threshold: total   CATE crosses +0.75
DATA <- Sys.getenv("AUTOC_DATA", unset = "data/derived/analysis_long.xlsx")

# ---- shared helpers (kept self-contained so this file runs on its own) -------
auc_trapezoid <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# fit LMM once, return BOTH total and chronic ATE (treatment - control), per year
fit_ate <- function(dataset, knot = KNOT) {
  dataset$time1 <- lspline(dataset$days, knot * 31)[, 1]
  dataset$time2 <- lspline(dataset$days, knot * 31)[, 2]
  fit <- lme(egfr ~ time1 * arm + time2 * arm, data = dataset,
             random = ~ 1 + time1 + time2 | id, control = lmeControl(opt = "optim"))
  fx <- fixef(fit); rf <- ranef(fit); nm <- names(fx)
  t1a <- nm[grepl("time1", nm) & grepl("arm", nm)]
  t2a <- nm[grepl("time2", nm) & grepl("arm", nm)]
  t1_no <- (fx["time1"] + rf[, "time1"]) * 365
  t1_ar <- (fx["time1"] + fx[t1a] + rf[, "time1"]) * 365
  t2_no <- (fx["time2"] + rf[, "time2"]) * 365
  t2_ar <- (fx["time2"] + fx[t2a] + rf[, "time2"]) * 365
  tot_no <- t1_no * (knot / 36) + t2_no * ((36 - knot) / 36)
  tot_ar <- t1_ar * (knot / 36) + t2_ar * ((36 - knot) / 36)
  c(total = mean(tot_ar) - mean(tot_no), chronic = mean(t2_ar) - mean(t2_no))
}

# UACR-ranked CATE curve: refit the LMM on the top-p% highest-UACR subgroup at
# each percentile cutoff, reading off total + chronic CATE.
uacr_curve <- function(df, knot = KNOT, n_pct = N_PCT) {
  base <- df[!duplicated(df$id), ]
  v <- sort(unique(base$uacr[!is.na(base$uacr)]))
  if (length(v) < n_pct) return(NULL)
  pct     <- quantile(v, probs = seq(0.01, 1, by = 0.01), na.rm = TRUE)
  cutoffs <- sort(pct)[1:n_pct]          # UACR: high value = high risk, bottom-up
  tot <- chr <- rep(NA_real_, n_pct)
  for (j in seq_along(cutoffs)) {
    sub <- df[df$uacr >= cutoffs[j], ]
    sub <- sub[!is.na(sub$egfr) & !is.na(sub$days) & !is.na(sub$arm), ]
    a <- tryCatch(fit_ate(sub, knot), error = function(e) c(total = NA, chronic = NA))
    tot[j] <- a["total"]; chr[j] <- a["chronic"]
  }
  data.frame(percentile = seq_len(n_pct), cov_value = as.numeric(cutoffs),
             total_cate = tot, chronic_cate = chr)
}

# first upward crossing of `target`, linearly interpolated in UACR units.
# Also flags a boundary crossing (target already met at the broadest subgroup
# or only at the very deepest one) and counts how many times the curve crosses
# target at all (a monotonicity / stability diagnostic).
cross_threshold <- function(curve, col, target) {
  y <- curve[[col]]; x <- curve$cov_value
  ncross <- sum(diff((y >= target) + 0L) != 0, na.rm = TRUE)
  idx <- which(y >= target)
  if (length(idx) == 0) return(c(thr = NA, crossed = 0, boundary = NA, ncross = ncross))
  j <- idx[1]
  boundary <- as.numeric(j == 1 || j == length(y))
  if (j == 1) return(c(thr = x[1], crossed = 1, boundary = 1, ncross = ncross))
  a0 <- y[j - 1]; a1 <- y[j]; c0 <- x[j - 1]; c1 <- x[j]
  thr <- if (is.na(a0) || is.na(a1) || a1 == a0) x[j]
         else c0 + (target - a0) / (a1 - a0) * (c1 - c0)
  c(thr = thr, crossed = 1, boundary = boundary, ncross = ncross)
}

# treatment effect above vs below a UACR threshold, for one slope.
eval_subgroups <- function(df, thr, slope) {
  ab <- df[df$uacr >= thr, ]; ab <- ab[!is.na(ab$egfr) & !is.na(ab$days) & !is.na(ab$arm), ]
  be <- df[df$uacr <  thr, ]; be <- be[!is.na(be$egfr) & !is.na(be$days) & !is.na(be$arm), ]
  aa <- tryCatch(fit_ate(ab), error = function(e) c(total = NA, chronic = NA))
  bb <- tryCatch(fit_ate(be), error = function(e) c(total = NA, chronic = NA))
  c(above = unname(aa[slope]), below = unname(bb[slope]),
    diff = unname(aa[slope] - bb[slope]),
    n_above = length(unique(ab$id)), n_below = length(unique(be$id)))
}

# ---- run --------------------------------------------------------------------
al <- read_excel(DATA)
al <- al[!is.na(al$egfr) & !is.na(al$days) & !is.na(al$arm) & !is.na(al$id) & !is.na(al$uacr), ]
cat("Pooled patients:", length(unique(al$id)), " rows:", nrow(al), "\n\n")

cu <- uacr_curve(al)
if (is.null(cu)) stop("Could not build UACR curve (too few distinct UACR values).")

autoc_chronic <- auc_trapezoid(N_PCT:1, cu$chronic_cate - cu$chronic_cate[1])
autoc_total   <- auc_trapezoid(N_PCT:1, cu$total_cate   - cu$total_cate[1])

tc <- cross_threshold(cu, "chronic_cate", TARGET_CHRONIC)
tt <- cross_threshold(cu, "total_cate",   TARGET_TOTAL)

sc <- if (tc["crossed"] == 1) {
  eval_subgroups(al, tc["thr"], "chronic")
} else c(above = NA, below = NA, diff = NA, n_above = NA, n_below = NA)
st <- if (tt["crossed"] == 1) {
  eval_subgroups(al, tt["thr"], "total")
} else c(above = NA, below = NA, diff = NA, n_above = NA, n_below = NA)

cat("================  PART 1: full-data point estimates  ================\n")
cat(sprintf("AUTOC(UACR) chronic slope : %.3f\n", autoc_chronic))
cat(sprintf("AUTOC(UACR) total   slope : %.3f\n\n", autoc_total))
cat(sprintf("CHRONIC threshold (CATE>%.2f): UACR = %s   [crossed=%d, boundary=%s, ncross=%d]\n",
            TARGET_CHRONIC, ifelse(is.na(tc["thr"]), "NA", sprintf("%.1f", tc["thr"])),
            tc["crossed"], tc["boundary"], tc["ncross"]))
cat(sprintf("  above (UACR>=thr): chronic ATE %.3f (n=%s) | below: %.3f (n=%s) | diff %.3f\n",
            sc["above"], sc["n_above"], sc["below"], sc["n_below"], sc["diff"]))
cat(sprintf("TOTAL   threshold (CATE>%.2f): UACR = %s   [crossed=%d, boundary=%s, ncross=%d]\n",
            TARGET_TOTAL, ifelse(is.na(tt["thr"]), "NA", sprintf("%.1f", tt["thr"])),
            tt["crossed"], tt["boundary"], tt["ncross"]))
cat(sprintf("  above (UACR>=thr): total ATE %.3f (n=%s) | below: %.3f (n=%s) | diff %.3f\n",
            st["above"], st["n_above"], st["below"], st["n_below"], st["diff"]))

# ---- save -------------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
point <- data.frame(
  slope       = c("chronic", "total"),
  target      = c(TARGET_CHRONIC, TARGET_TOTAL),
  autoc_uacr  = c(autoc_chronic, autoc_total),
  threshold   = c(unname(tc["thr"]),      unname(tt["thr"])),
  crossed     = c(unname(tc["crossed"]),  unname(tt["crossed"])),
  boundary    = c(unname(tc["boundary"]), unname(tt["boundary"])),
  ncross      = c(unname(tc["ncross"]),   unname(tt["ncross"])),
  ate_above   = c(unname(sc["above"]),    unname(st["above"])),
  ate_below   = c(unname(sc["below"]),    unname(st["below"])),
  ate_diff    = c(unname(sc["diff"]),     unname(st["diff"])),
  n_above     = c(unname(sc["n_above"]),  unname(st["n_above"])),
  n_below     = c(unname(sc["n_below"]),  unname(st["n_below"])),
  row.names = NULL)
write.csv(point,          "results/cv_point.csv",       row.names = FALSE)
write.csv(cu,             "results/cv_point_curve.csv",  row.names = FALSE)
cat("\nWrote results/cv_point.csv and results/cv_point_curve.csv\n")
