# =============================================================================
# 06b_combine_bootcv.R  --  Variance-component step of bootstrap-CV
# (Cai, ..., Tian 2025, arXiv:2307.00260, Algorithm 2 step 11).
#
# Reads frag_bootcv_*.csv and fits, for each estimand, the random-effects model
#
#     theta*_bk = theta0 + eps*_b + eps_bk ,   eps*_b ~ (0, sigma^2_BT)
#                                              eps_bk ~ (0, tau^2_0)
#
# via the paper's closed-form moment estimators:
#
#   tau^2_0    = sum_b sum_k (theta*_bk - mean_b)^2 / ((B_CV - 1) * B_BOOT)
#   sigma^2_BT = sum_b (mean_b - grand_mean)^2 / (B_BOOT - 1)  -  tau^2_0 / B_CV
#
# sigma^2_BT is the BETWEEN-bootstrap variance = the sampling uncertainty we
# want. tau^2_0 is the WITHIN-bootstrap variance = noise from random train/test
# splitting, a nuisance. Subtracting tau^2_0/B_CV removes the splitting noise
# that contaminates the between-bootstrap spread -- this is what makes a small
# B_CV (10-20) sufficient instead of 200+.
#
#   95% CI  =  Err^CV_m  +/-  1.96 * sigma_BT
#
# The POINT estimate Err^CV_m comes from the unweighted repeated-split CV at the
# NOMINAL training fraction (05_autoc_cv.R); this script supplies only the
# standard error. Pass it via AUTOC_MCCV. (m_adj is used inside the bootstrap
# purely to keep the *effective* training size near m -- it is not the estimand.)
#
# Remark 3 deflation: bootstrapped training sets hold only 0.632*m_adj distinct
# patients, inflating the variance, so we also report
#     sigma^2_adj = sigma^2_BT * (1 - 0.368 * m_adj / n).
#
# Env vars:
#   AUTOC_OUT     where frag_bootcv_*.csv are  (default $SCRATCH, else ".")
#   AUTOC_FINAL   where to write finals        (default ~)
#   AUTOC_MCCV    path to autoc_cv_raw.csv from 05b (for the point estimate)
#   N_PATIENTS    n, for the Remark-3 deflation (default: read from fragments)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(writexl) })

OUT   <- Sys.getenv("AUTOC_OUT",   unset = Sys.getenv("SCRATCH", unset = "."))
FINAL <- Sys.getenv("AUTOC_FINAL", unset = Sys.getenv("HOME", unset = "."))
MCCV  <- Sys.getenv("AUTOC_MCCV",  unset = "")

files <- list.files(OUT, pattern = "^frag_bootcv_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No frag_bootcv_*.csv fragments in ", OUT)
bc <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
cat("Combined", nrow(bc), "(bootstrap x CV) units from", length(files), "fragments.\n")
cat("Distinct bootstraps:", n_distinct(bc$boot),
    "| mean CV splits per bootstrap:", round(nrow(bc) / n_distinct(bc$boot), 1), "\n")

M_ADJ <- bc$m_adj[1]
N_PT  <- as.numeric(Sys.getenv("N_PATIENTS", unset = NA))

# ---- the variance-component estimator ---------------------------------------
# Uses only (b, k) units where the estimand is non-missing (threshold found and
# both subgroup fits converged). Bootstraps contributing <2 usable splits cannot
# inform tau^2 and are dropped from that term.
bootcv_se <- function(df, col) {
  d <- df[!is.na(df[[col]]), c("boot", col)]
  names(d) <- c("boot", "y")
  if (nrow(d) < 10) return(NULL)

  per_b <- d %>% group_by(boot) %>%
    summarise(nb = n(), mb = mean(y), ss = sum((y - mean(y))^2), .groups = "drop")

  B_BOOT   <- nrow(per_b)
  B_CV_eff <- mean(per_b$nb)                       # harmonised inner count

  # within-bootstrap (splitting) variance
  denom  <- sum(per_b$nb - 1)
  tau2   <- if (denom > 0) sum(per_b$ss) / denom else NA_real_

  # between-bootstrap spread, then subtract the splitting-noise contribution
  grand  <- mean(per_b$mb)
  s2_raw <- sum((per_b$mb - grand)^2) / (B_BOOT - 1)
  sigma2 <- s2_raw - tau2 / B_CV_eff

  truncated <- sigma2 <= 0
  if (truncated) sigma2 <- s2_raw / B_CV_eff       # conservative floor

  # Remark 3: deflate for the 0.632 distinct-patient loss in bootstrapped train
  defl <- if (!is.na(N_PT)) (1 - 0.368 * M_ADJ / N_PT) else NA_real_
  sigma2_adj <- if (!is.na(defl)) sigma2 * defl else NA_real_

  data.frame(
    estimand   = col,
    boot_mean  = grand,
    B_BOOT     = B_BOOT,
    B_CV_eff   = B_CV_eff,
    tau2_within  = tau2,
    s2_between_raw = s2_raw,
    sigma2_BT  = sigma2,
    se_BT      = sqrt(sigma2),
    sigma2_adj = sigma2_adj,
    se_adj     = if (!is.na(sigma2_adj)) sqrt(sigma2_adj) else NA_real_,
    # Remark 2: optimal B_CV ~ tau2/sigma2_BT. If this exceeds the B_CV you ran,
    # the SE is noisier than it needs to be -- raise B_CV next time.
    optimal_B_CV = tau2 / sigma2,
    var_floor_hit = truncated,
    n_units    = nrow(d),
    row.names  = NULL)
}

ESTIMANDS <- c("chronic_above", "chronic_below", "chronic_diff",
               "total_above",   "total_below",   "total_diff")
var_tab <- do.call(rbind, lapply(ESTIMANDS, function(e) bootcv_se(bc, e)))

# ---- point estimate from the unweighted MC-CV (05b) -------------------------
point <- setNames(rep(NA_real_, length(ESTIMANDS)), ESTIMANDS)
point_src <- "bootstrap mean (FALLBACK - supply AUTOC_MCCV for the correct Err^CV_m)"
if (nzchar(MCCV) && file.exists(MCCV)) {
  mc <- read.csv(MCCV, stringsAsFactors = FALSE)
  map <- c(chronic_above = "te_chronic_above", chronic_below = "te_chronic_below",
           chronic_diff  = "te_chronic_diff",  total_above   = "te_total_above",
           total_below   = "te_total_below",   total_diff    = "te_total_diff")
  for (e in ESTIMANDS) if (!is.null(mc[[map[[e]]]]))
    point[[e]] <- mean(mc[[map[[e]]]], na.rm = TRUE)
  point_src <- paste("MC-CV mean from", MCCV)
  cat("Point estimate source:", point_src, "\n")
} else {
  for (e in ESTIMANDS) point[[e]] <- var_tab$boot_mean[var_tab$estimand == e]
  cat("WARNING: AUTOC_MCCV not supplied; using bootstrap mean as point estimate.\n")
}

# ---- final CIs --------------------------------------------------------------
final <- var_tab %>%
  mutate(point_estimate = as.numeric(point[estimand]),
         lower_95 = point_estimate - 1.96 * se_BT,
         upper_95 = point_estimate + 1.96 * se_BT,
         lower_95_adj = point_estimate - 1.96 * se_adj,
         upper_95_adj = point_estimate + 1.96 * se_adj) %>%
  select(estimand, point_estimate, se_BT, lower_95, upper_95,
         se_adj, lower_95_adj, upper_95_adj, everything())

# ---- side-by-side with the naive MC-CV percentile interval ------------------
comparison <- NULL
if (nzchar(MCCV) && file.exists(MCCV)) {
  mc <- read.csv(MCCV, stringsAsFactors = FALSE)
  map <- c(chronic_above = "te_chronic_above", chronic_below = "te_chronic_below",
           chronic_diff  = "te_chronic_diff",  total_above   = "te_total_above",
           total_below   = "te_total_below",   total_diff    = "te_total_diff")
  comparison <- do.call(rbind, lapply(ESTIMANDS, function(e) {
    x <- mc[[map[[e]]]]
    if (is.null(x)) return(NULL)
    naive_lo <- quantile(x, 0.025, na.rm = TRUE, names = FALSE)
    naive_hi <- quantile(x, 0.975, na.rm = TRUE, names = FALSE)
    f <- final[final$estimand == e, ]
    data.frame(estimand = e, point = mean(x, na.rm = TRUE),
               naive_mccv_lower = naive_lo, naive_mccv_upper = naive_hi,
               naive_width = naive_hi - naive_lo,
               bootcv_lower = f$lower_95, bootcv_upper = f$upper_95,
               bootcv_width = f$upper_95 - f$lower_95,
               width_ratio = (f$upper_95 - f$lower_95) / (naive_hi - naive_lo),
               row.names = NULL)
  }))
}

# ---- console ----------------------------------------------------------------
cat("\n=========  Bootstrap-CV variance components  =========\n")
print(var_tab[, c("estimand", "tau2_within", "sigma2_BT", "se_BT", "optimal_B_CV", "n_units")],
      row.names = FALSE, digits = 4)
cat("\n=========  HEADLINE: subgroup CATEs with valid 95% CI  =========\n")
print(final[, c("estimand", "point_estimate", "se_BT", "lower_95", "upper_95")],
      row.names = FALSE, digits = 4)
if (!is.null(comparison)) {
  cat("\n=========  Naive MC-CV percentile CI vs bootstrap-CV CI  =========\n")
  print(comparison[, c("estimand", "naive_width", "bootcv_width", "width_ratio")],
        row.names = FALSE, digits = 4)
  cat("\n(width_ratio >> 1 means the naive percentile interval was understating\n",
      " uncertainty, as expected -- it measures splitting noise, not sampling error.)\n")
}

# ---- write ------------------------------------------------------------------
sheets <- list(headline_ci      = final,
               variance_comps   = var_tab,
               raw_bootcv_units = bc)
if (!is.null(comparison)) sheets$mccv_vs_bootcv <- comparison

res_file <- file.path(FINAL, "autoc_bootcv_results.xlsx")
write_xlsx(sheets, res_file)
write.csv(bc, file.path(FINAL, "autoc_bootcv_raw.csv"), row.names = FALSE)
cat("\nWrote:\n ", res_file, "\n ", file.path(FINAL, "autoc_bootcv_raw.csv"), "\n")
cat("Point estimate source:", point_src, "\n")
