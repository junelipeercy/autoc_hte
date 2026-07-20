# =============================================================================
# 03b_combine_bootstrap.R  --  Stitch one cohort's array fragments into results.
#
# Run ONCE per cohort after its array finishes:
#     COHORT=SPRINT Rscript 03b_combine_bootstrap.R
# (light -- fine on an sh_dev node; not on the login node for big jobs.)
#
# Reads frag_{scalars,curves}_<COHORT>_*.csv from the fragment dir and writes:
#   autoc_bootstrap_results_<COHORT>.xlsx   (summaries + merged raw)
#   autoc_curves_for_plot_<COHORT>.csv      (mean + 95% CI ribbon)
#
# Env vars:
#   COHORT       POOLED | SPRINT | ACCORD | AASK | MDRD   (default POOLED)
#   AUTOC_OUT    where the fragments are (default $SCRATCH, else ".")
#   AUTOC_FINAL  where to write finals   (default ~)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(writexl) })

COHORT <- toupper(Sys.getenv("COHORT", "POOLED"))
OUT    <- Sys.getenv("AUTOC_OUT",   unset = Sys.getenv("SCRATCH", unset = "."))
FINAL  <- Sys.getenv("AUTOC_FINAL", unset = Sys.getenv("HOME", unset = "."))

read_all <- function(kind) {
  pat <- sprintf("^frag_%s_%s_.*\\.csv$", kind, COHORT)
  files <- list.files(OUT, pattern = pat, full.names = TRUE)
  if (length(files) == 0) stop("No fragments matching ", pat, " in ", OUT)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}
scalars <- read_all("scalars")
curves  <- read_all("curves")
cat(COHORT, ": combined", nrow(scalars), "bootstraps,", nrow(curves), "curve rows.\n")

ci <- function(x) c(mean = mean(x, na.rm = TRUE),
                    lower_95 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
                    upper_95 = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
                    n_ok = sum(!is.na(x)))

# 1. AUTOC scores per biomarker x slope
autoc_summary <- do.call(rbind, lapply(
  c("autoc_chronic_egfr","autoc_chronic_uacr","autoc_chronic_kfre",
    "autoc_total_egfr","autoc_total_uacr","autoc_total_kfre"),
  function(col) data.frame(cohort = COHORT, metric = col, t(ci(scalars[[col]])))))

# 2. Pairwise AUTOC differences
pw <- function(a, b, lab) data.frame(cohort = COHORT, comparison = lab, t(ci(scalars[[a]] - scalars[[b]])))
pairwise_summary <- rbind(
  pw("autoc_chronic_uacr","autoc_chronic_kfre","chronic: UACR-KFRE"),
  pw("autoc_chronic_uacr","autoc_chronic_egfr","chronic: UACR-eGFR"),
  pw("autoc_chronic_kfre","autoc_chronic_egfr","chronic: KFRE-eGFR"),
  pw("autoc_total_uacr","autoc_total_kfre","total: UACR-KFRE"),
  pw("autoc_total_uacr","autoc_total_egfr","total: UACR-eGFR"),
  pw("autoc_total_kfre","autoc_total_egfr","total: KFRE-eGFR"))

# 3. Winner tally
winner_summary <- rbind(
  data.frame(cohort = COHORT, slope = "chronic", as.data.frame(t(prop.table(table(
    factor(scalars$winner_chronic, c("eGFR","UACR","KFRE"))))))),
  data.frame(cohort = COHORT, slope = "total", as.data.frame(t(prop.table(table(
    factor(scalars$winner_total, c("eGFR","UACR","KFRE"))))))))

# 4. UACR thresholds + crossing rate
threshold_summary <- data.frame(
  cohort = COHORT,
  slope = c("chronic (CATE>1.0)", "total (CATE>0.75)"),
  crossing_rate = c(mean(scalars$crossed_chronic), mean(scalars$crossed_total)),
  rbind(ci(scalars$thr_chronic[scalars$crossed_chronic == 1]),
        ci(scalars$thr_total[scalars$crossed_total == 1])))

# 5. Subgroup CATEs (above/below) + difference (HTE test) + subgroup sizes
subgroup_summary <- do.call(rbind, lapply(
  c("sg_chronic_above","sg_chronic_below","sg_chronic_diff",
    "n_chronic_above","n_chronic_below",
    "sg_total_above","sg_total_below","sg_total_diff",
    "n_total_above","n_total_below"),
  function(col) data.frame(cohort = COHORT, metric = col, t(ci(scalars[[col]])))))

# 6. Overall pooled ATE (headline)
overall_summary <- rbind(
  data.frame(cohort = COHORT, metric = "overall_chronic", t(ci(scalars$overall_chronic))),
  data.frame(cohort = COHORT, metric = "overall_total",   t(ci(scalars$overall_total))))

# 7. AUTOC curve ribbon
curve_plot <- curves %>%
  group_by(boot, biomarker) %>%
  mutate(chronic_diff = chronic_cate - chronic_cate[percentile == 1],
         total_diff   = total_cate   - total_cate[percentile == 1]) %>%
  ungroup() %>%
  group_by(biomarker, percentile) %>%
  summarise(cohort         = COHORT,
            cov_value      = mean(cov_value, na.rm = TRUE),
            chronic_mean   = mean(chronic_diff, na.rm = TRUE),
            chronic_lower  = quantile(chronic_diff, 0.025, na.rm = TRUE),
            chronic_upper  = quantile(chronic_diff, 0.975, na.rm = TRUE),
            total_mean     = mean(total_diff, na.rm = TRUE),
            total_lower    = quantile(total_diff, 0.025, na.rm = TRUE),
            total_upper    = quantile(total_diff, 0.975, na.rm = TRUE),
            .groups = "drop")

# ---- write finals -----------------------------------------------------------
res_file   <- file.path(FINAL, paste0("autoc_bootstrap_results_", COHORT, ".xlsx"))
curve_file <- file.path(FINAL, paste0("autoc_curves_for_plot_", COHORT, ".csv"))
write_xlsx(list(autoc_scores     = autoc_summary,
                autoc_pairwise   = pairwise_summary,
                winner_tally     = winner_summary,
                uacr_thresholds  = threshold_summary,
                subgroup_cates   = subgroup_summary,
                overall_ate      = overall_summary,
                raw_bootstraps   = scalars),
           res_file)
write.csv(curve_plot, curve_file, row.names = FALSE)

cat("\nWrote:\n ", res_file, "\n ", curve_file, "\n")
cat("\n", COHORT, "UACR crossing rates -> chronic:",
    round(mean(scalars$crossed_chronic), 3), " total:",
    round(mean(scalars$crossed_total), 3), "\n")
