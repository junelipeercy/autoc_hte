# =============================================================================
# 05b_combine_cv.R  --  PART 3: stitch the CV fragments into results + CIs.
#
# Run ONCE after the 05_autoc_cv.R array finishes:
#     Rscript 05b_combine_cv.R
# (light -- fine on an sh_dev node.)
#
# Reads frag_cv_*.csv from the fragment dir and writes:
#   autoc_cv_results.xlsx   (all summary sheets + merged raw iterations)
#   autoc_cv_raw.csv        (one row per CV iteration, for custom plots)
#
# Headline (item 3): held-out TEST treatment effect above vs below the UACR
# threshold, with 95% CIs, for chronic and total slope.
# Also reports the crossing/success rate, threshold distribution, and the
# optimism (TRAIN in-sample minus TEST held-out) that quantifies AUTOC's bias.
#
# Env vars:
#   AUTOC_OUT     where fragments are        (default $SCRATCH, else ".")
#   AUTOC_FINAL   where to write finals      (default ~)
#   AUTOC_POINT   optional path to results/cv_point.csv from 04_cv_point.R;
#                 if present, its full-data numbers are added as the optimistic
#                 tier of the bias comparison.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(writexl) })

OUT   <- Sys.getenv("AUTOC_OUT",   unset = Sys.getenv("SCRATCH", unset = "."))
FINAL <- Sys.getenv("AUTOC_FINAL", unset = Sys.getenv("HOME", unset = "."))
POINT <- Sys.getenv("AUTOC_POINT", unset = "")

files <- list.files(OUT, pattern = "^frag_cv_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No frag_cv_*.csv fragments in ", OUT)
cv <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
cat("Combined", nrow(cv), "CV iterations from", length(files), "fragments.\n")

# percentile CI over CV iterations (NA-robust)
ci <- function(x) c(mean = mean(x, na.rm = TRUE),
                    median = median(x, na.rm = TRUE),
                    lower_95 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
                    upper_95 = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
                    n_ok = sum(!is.na(x)))
row_ci <- function(label, x) data.frame(metric = label, t(ci(x)), row.names = NULL)

# "usable" iterations per slope: threshold crossed AND both TEST subgroup fits
# converged. These are the ones that feed the effect summaries.
ok_chronic <- cv$crossed_chronic == 1 & cv$ok_te_chronic_above == 1 & cv$ok_te_chronic_below == 1
ok_total   <- cv$crossed_total   == 1 & cv$ok_te_total_above   == 1 & cv$ok_te_total_below   == 1
ok_chronic[is.na(ok_chronic)] <- FALSE
ok_total[is.na(ok_total)]     <- FALSE

# ---- 1. success / crossing / usability rates (a PRIMARY result) -------------
rate_summary <- data.frame(
  slope          = c("chronic (CATE>1.00)", "total (CATE>0.75)"),
  n_iterations   = c(nrow(cv), nrow(cv)),
  crossing_rate  = c(mean(cv$crossed_chronic, na.rm = TRUE), mean(cv$crossed_total, na.rm = TRUE)),
  interior_rate  = c(mean(cv$crossed_chronic == 1 & cv$boundary_chronic == 0, na.rm = TRUE),
                     mean(cv$crossed_total   == 1 & cv$boundary_total   == 0, na.rm = TRUE)),
  usable_rate    = c(mean(ok_chronic), mean(ok_total)),
  n_usable       = c(sum(ok_chronic), sum(ok_total)),
  mean_ncross    = c(mean(cv$ncross_chronic, na.rm = TRUE), mean(cv$ncross_total, na.rm = TRUE)),
  row.names = NULL)

# ---- 2. threshold distribution (stability of the discovered cutoff) ---------
threshold_summary <- rbind(
  row_ci("thr_chronic (UACR)", cv$thr_chronic[ok_chronic]),
  row_ci("thr_total (UACR)",   cv$thr_total[ok_total]))

# ---- 3. HEADLINE: held-out TEST treatment effect above/below threshold ------
test_summary <- rbind(
  row_ci("chronic: TEST ATE above", cv$te_chronic_above[ok_chronic]),
  row_ci("chronic: TEST ATE below", cv$te_chronic_below[ok_chronic]),
  row_ci("chronic: TEST diff (above-below)", cv$te_chronic_diff[ok_chronic]),
  row_ci("total: TEST ATE above",   cv$te_total_above[ok_total]),
  row_ci("total: TEST ATE below",   cv$te_total_below[ok_total]),
  row_ci("total: TEST diff (above-below)",   cv$te_total_diff[ok_total]))

# ---- 4. in-sample TRAIN effect at the same threshold (optimistic tier) ------
train_summary <- rbind(
  row_ci("chronic: TRAIN ATE above", cv$tr_chronic_above[ok_chronic]),
  row_ci("chronic: TRAIN ATE below", cv$tr_chronic_below[ok_chronic]),
  row_ci("chronic: TRAIN diff (above-below)", cv$tr_chronic_diff[ok_chronic]),
  row_ci("total: TRAIN ATE above",   cv$tr_total_above[ok_total]),
  row_ci("total: TRAIN ATE below",   cv$tr_total_below[ok_total]),
  row_ci("total: TRAIN diff (above-below)",   cv$tr_total_diff[ok_total]))

# ---- 5. OPTIMISM / BIAS: paired TRAIN - TEST within each iteration -----------
# Positive optimism = the in-sample effect overstates the honest held-out effect,
# i.e. the AUTOC discovery process is optimistically biased by that much.
opt_chronic <- cv$tr_chronic_diff[ok_chronic] - cv$te_chronic_diff[ok_chronic]
opt_total   <- cv$tr_total_diff[ok_total]     - cv$te_total_diff[ok_total]
optimism_summary <- rbind(
  row_ci("chronic: optimism (TRAINdiff - TESTdiff)", opt_chronic),
  row_ci("total: optimism (TRAINdiff - TESTdiff)",   opt_total))

# ---- 6. subgroup sizes in the held-out fold (stability check) ---------------
size_summary <- rbind(
  row_ci("chronic: n TEST above", cv$n_te_chronic_above[ok_chronic]),
  row_ci("chronic: n TEST below", cv$n_te_chronic_below[ok_chronic]),
  row_ci("total: n TEST above",   cv$n_te_total_above[ok_total]),
  row_ci("total: n TEST below",   cv$n_te_total_below[ok_total]))

# ---- optional: fold in the full-data point estimate (most optimistic tier) --
point_summary <- NULL
if (nzchar(POINT) && file.exists(POINT)) {
  pt <- read.csv(POINT, stringsAsFactors = FALSE)
  point_summary <- data.frame(
    slope = pt$slope, threshold = pt$threshold,
    fulldata_ate_above = pt$ate_above, fulldata_ate_below = pt$ate_below,
    fulldata_diff = pt$ate_diff, row.names = NULL)
  cat("Included full-data point estimates from", POINT, "\n")
}

# ---- console headline -------------------------------------------------------
cat("\n=========  CV usability  =========\n"); print(rate_summary, row.names = FALSE)
cat("\n=========  Held-out TEST effects (headline)  =========\n"); print(test_summary, row.names = FALSE)
cat("\n=========  Optimism (bias of AUTOC discovery)  =========\n"); print(optimism_summary, row.names = FALSE)

# ---- write finals -----------------------------------------------------------
sheets <- list(
  usability_rates  = rate_summary,
  thresholds       = threshold_summary,
  test_effects     = test_summary,     # <- item 3 headline
  train_effects    = train_summary,
  optimism_bias    = optimism_summary,
  subgroup_sizes   = size_summary,
  raw_iterations   = cv)
if (!is.null(point_summary)) sheets$fulldata_point <- point_summary

res_file <- file.path(FINAL, "autoc_cv_results.xlsx")
raw_file <- file.path(FINAL, "autoc_cv_raw.csv")
write_xlsx(sheets, res_file)
write.csv(cv, raw_file, row.names = FALSE)
cat("\nWrote:\n ", res_file, "\n ", raw_file, "\n")
