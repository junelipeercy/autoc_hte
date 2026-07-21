# =============================================================================
# 07_cv_figures.R  --  Dissertation figures + assembled tables for the UACR
# threshold cross-validation (steps 04-06b).
#
# Reads (from results/):
#   cv_point.csv               full-data point thresholds (step 04)
#   autoc_cv_raw.csv           1000 MC-CV iterations (step 05/05b)
#   autoc_cv_results.xlsx       MC-CV summaries (step 05b)
#   autoc_bootcv_results.xlsx   bootstrap-CV valid CIs (step 06b)
#
# Writes (to figures/ and results/):
#   fig1_forest_subgroup_cate.(pdf|png)   headline: above/below CATE + valid CI
#   fig2_threshold_stability.(pdf|png)    discovered UACR threshold across runs
#   fig3_variance_decomposition.(pdf|png) splitting noise vs sampling error
#   dissertation_tables.xlsx              all tables, one workbook, print-ready
#
# Design: diverging colour by polarity (protective = blue, harmful = orange),
# CVD-safe hues, thin marks, direct labels, neutral grid. Run locally.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(writexl); library(dplyr); library(tidyr); library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)
UNIT <- expression(paste("Chronic / total eGFR slope ATE  (mL/min/1.73",m^2,"/yr)"))

# wrap long subtitle/caption strings so they don't clip at the panel edge
wrap <- function(s, n = 95) paste(strwrap(s, width = n), collapse = "\n")

# CVD-safe palette (dataviz reference instance)
COL_ABOVE <- "#2a78d6"   # blue   -- above threshold (protective)
COL_BELOW <- "#eb6834"   # orange -- below threshold
COL_DIFF  <- "#4a3aa7"   # violet -- above-below difference
COL_POINT <- "#0b0b0b"
COL_MUTED <- "#8a8987"
COL_GRID  <- "#e6e5e2"

theme_diss <- function(base = 12) {
  theme_minimal(base_size = base, base_family = "") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
          axis.title = element_text(colour = "#2b2b2b"),
          axis.text  = element_text(colour = "#2b2b2b"),
          strip.text = element_text(face = "bold", colour = "#2b2b2b", size = base),
          plot.title = element_text(face = "bold", size = base + 3),
          plot.subtitle = element_text(colour = "#52514e", size = base - 1),
          plot.caption  = element_text(colour = COL_MUTED, size = base - 3, hjust = 0),
          legend.position = "none",
          plot.margin = margin(14, 18, 12, 14))
}

# ---- load -------------------------------------------------------------------
point <- read.csv("results/cv_point.csv", stringsAsFactors = FALSE)
raw   <- read.csv("results/autoc_cv_raw.csv", stringsAsFactors = FALSE)
h     <- as.data.frame(read_excel("results/autoc_bootcv_results.xlsx", "headline_ci"))
vcomp <- as.data.frame(read_excel("results/autoc_bootcv_results.xlsx", "variance_comps"))
cmp   <- as.data.frame(read_excel("results/autoc_bootcv_results.xlsx", "mccv_vs_bootcv"))

# =============================================================================
# FIGURE 1 -- Forest plot of held-out subgroup CATEs with VALID boot-CV 95% CI,
# with the (wider, invalid) naive percentile interval drawn faintly behind to
# make the "naive overstates uncertainty" point visible.
# =============================================================================
lab_map <- c(
  chronic_above = "Above threshold",  chronic_below = "Below threshold",  chronic_diff = "Difference (above - below)",
  total_above   = "Above threshold",  total_below   = "Below threshold",  total_diff   = "Difference (above - below)")
slope_of <- function(e) ifelse(grepl("^chronic", e), "Chronic slope  (threshold: CATE > 1.0)",
                                                      "Total slope  (threshold: CATE > 0.75)")
role_of  <- function(e) sub(".*_", "", e)

fig1_df <- h %>%
  transmute(estimand,
            slope = slope_of(estimand), role = role_of(estimand),
            label = lab_map[estimand],
            point = point_estimate,
            lo = lower_95_adj, hi = upper_95_adj) %>%          # Remark-3 adjusted CI
  left_join(cmp %>% transmute(estimand, naive_lo = naive_mccv_lower, naive_hi = naive_mccv_upper),
            by = "estimand") %>%
  mutate(label = factor(label, levels = c("Difference (above - below)", "Below threshold", "Above threshold")),
         col = c(above = COL_ABOVE, below = COL_BELOW, diff = COL_DIFF)[role])

fig1 <- ggplot(fig1_df, aes(y = label)) +
  geom_vline(xintercept = 0, colour = "#b0afab", linewidth = 0.4) +
  geom_linerange(aes(xmin = naive_lo, xmax = naive_hi), colour = COL_MUTED,
                 linewidth = 3.2, alpha = 0.28) +                               # naive (invalid)
  geom_linerange(aes(xmin = lo, xmax = hi, colour = col), linewidth = 1.4) +    # valid CI
  geom_point(aes(x = point, colour = col), size = 3.1) +
  geom_text(aes(x = point, label = sprintf("%.2f", point)), vjust = -1.1, size = 3.5, colour = COL_POINT) +
  geom_text(aes(x = hi, label = sprintf("[%.2f, %.2f]", lo, hi)), hjust = -0.12,
            size = 2.9, colour = "#52514e") +
  facet_wrap(~ slope, ncol = 1, scales = "free_y") +
  scale_colour_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.20))) +
  labs(title = "Held-out treatment effect above vs below the UACR threshold",
       subtitle = wrap("Point = cross-validated CATE; thick bar = valid bootstrap-CV 95% CI; faint grey bar = naive percentile interval"),
       x = UNIT, y = NULL,
       caption = wrap("1000-iteration Monte-Carlo CV for the point estimate; bootstrap-CV (Cai et al. 2025) for the CI. Positive = intensive BP control preserves eGFR.", 110)) +
  theme_diss()
ggsave("figures/fig1_forest_subgroup_cate.pdf", fig1, width = 9.0, height = 6.4)
ggsave("figures/fig1_forest_subgroup_cate.png", fig1, width = 9.0, height = 6.4, dpi = 300)

# =============================================================================
# FIGURE 2 -- Stability of the DISCOVERED threshold across the 1000 CV runs.
# Shows that the method reliably finds *a* cutpoint but its exact value is
# loosely determined (esp. total slope).
# =============================================================================
thr_df <- bind_rows(
  data.frame(slope = "Chronic slope", thr = raw$thr_chronic[raw$crossed_chronic == 1]),
  data.frame(slope = "Total slope",   thr = raw$thr_total[raw$crossed_total == 1])) %>%
  filter(!is.na(thr))

thr_stats <- thr_df %>% group_by(slope) %>%
  summarise(med = median(thr),
            lo = quantile(thr, 0.025), hi = quantile(thr, 0.975), .groups = "drop")
pt_lines <- data.frame(slope = c("Chronic slope", "Total slope"),
                       full = point$threshold[match(c("chronic","total"), point$slope)])
col_slope <- c("Chronic slope" = COL_ABOVE, "Total slope" = COL_BELOW)

fig2 <- ggplot(thr_df, aes(thr, fill = slope)) +
  geom_histogram(bins = 45, colour = "white", linewidth = 0.15) +
  geom_vline(data = thr_stats, aes(xintercept = med), colour = COL_POINT, linewidth = 0.6) +
  geom_vline(data = thr_stats, aes(xintercept = lo), colour = COL_MUTED, linetype = "22", linewidth = 0.4) +
  geom_vline(data = thr_stats, aes(xintercept = hi), colour = COL_MUTED, linetype = "22", linewidth = 0.4) +
  geom_vline(data = pt_lines, aes(xintercept = full), colour = COL_DIFF, linetype = "42", linewidth = 0.6) +
  geom_text(data = thr_stats, aes(x = med, y = Inf, label = sprintf("median %.0f", med)),
            vjust = 1.6, hjust = -0.08, size = 3, colour = COL_POINT) +
  facet_wrap(~ slope, ncol = 1, scales = "free") +
  scale_fill_manual(values = col_slope) +
  scale_x_continuous(limits = c(0, 300), oob = scales::oob_squish) +
  labs(title = "Where does the UACR threshold land across 1000 CV splits?",
       subtitle = "Solid line = CV median; dashed grey = 95% range; dotted violet = full-data point estimate",
       x = "Discovered UACR threshold  (mg/g)", y = "CV iterations",
       caption = wrap("Truncated at 300 mg/g for display. A reliable separating cutpoint is always found, but its exact value is loosely determined (esp. total slope).", 110)) +
  theme_diss()
ggsave("figures/fig2_threshold_stability.pdf", fig2, width = 8.6, height = 6.0)
ggsave("figures/fig2_threshold_stability.png", fig2, width = 8.6, height = 6.0, dpi = 300)

# =============================================================================
# FIGURE 3 -- Variance decomposition: splitting noise (tau^2) vs sampling error
# (sigma^2_BT). The methods point -- why the naive percentile interval, which
# tracks tau, overstates the uncertainty of the averaged CV estimate (sigma_BT).
# =============================================================================
vc_df <- vcomp %>%
  transmute(estimand,
            slope = slope_of(estimand), label = lab_map[estimand],
            `Splitting noise (SD, within)` = sqrt(tau2_within),
            `Sampling error (SE, between)` = se_BT) %>%
  filter(role_of(estimand) != "diff") %>%
  pivot_longer(c(`Splitting noise (SD, within)`, `Sampling error (SE, between)`),
               names_to = "component", values_to = "sd") %>%
  mutate(label = factor(label, levels = c("Below threshold", "Above threshold")),
         component = factor(component, levels = c("Splitting noise (SD, within)", "Sampling error (SE, between)")))

fig3 <- ggplot(vc_df, aes(sd, label, fill = component)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", sd)), position = position_dodge(width = 0.7),
            hjust = -0.2, size = 3, colour = "#2b2b2b") +
  facet_wrap(~ slope, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Splitting noise (SD, within)" = COL_MUTED,
                               "Sampling error (SE, between)" = COL_ABOVE)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Splitting noise dominates sampling error",
       subtitle = "The naive percentile interval tracks the grey bar; the valid CI tracks the blue bar",
       x = "Standard deviation of the subgroup CATE  (mL/min/1.73m2/yr)", y = NULL,
       caption = wrap("Splitting noise (within-bootstrap, tau) = spread of a single-split estimate. Sampling error (between-bootstrap, sigma_BT) = uncertainty of the averaged CV estimate.", 110)) +
  theme_diss() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 10))
ggsave("figures/fig3_variance_decomposition.pdf", fig3, width = 8.6, height = 5.4)
ggsave("figures/fig3_variance_decomposition.png", fig3, width = 8.6, height = 5.4, dpi = 300)

# =============================================================================
# TABLES -- assembled into one print-ready workbook.
# =============================================================================
rnd <- function(x, d = 2) round(as.numeric(x), d)

# Table 1: full-data point thresholds + in-sample subgroup effects (step 04)
t1 <- point %>% transmute(
  Slope = tools::toTitleCase(slope),
  `Target CATE` = target,
  `AUTOC (UACR)` = rnd(autoc_uacr, 1),
  `UACR threshold (mg/g)` = rnd(threshold, 1),
  `ATE above` = rnd(ate_above), `ATE below` = rnd(ate_below),
  `n above` = n_above, `n below` = n_below)

# Table 2: CV reliability -- success rate, threshold stability, optimism
usab <- as.data.frame(read_excel("results/autoc_cv_results.xlsx", "usability_rates"))
opt  <- as.data.frame(read_excel("results/autoc_cv_results.xlsx", "optimism_bias"))
t2 <- data.frame(
  Slope = c("Chronic", "Total"),
  `Success rate` = rnd(usab$crossing_rate, 3),
  `Usable rate`  = rnd(usab$usable_rate, 3),
  `Threshold median (mg/g)` = rnd(thr_stats$med[match(c("Chronic slope","Total slope"), thr_stats$slope)], 0),
  `Threshold 95% low`  = rnd(thr_stats$lo[match(c("Chronic slope","Total slope"), thr_stats$slope)], 0),
  `Threshold 95% high` = rnd(thr_stats$hi[match(c("Chronic slope","Total slope"), thr_stats$slope)], 0),
  `Optimism (train-test)` = rnd(opt$mean, 3),
  check.names = FALSE)

# Table 3: HEADLINE -- held-out subgroup CATEs with valid boot-CV 95% CI
t3 <- h %>% transmute(
  Slope = slope_of(estimand),
  Subgroup = lab_map[estimand],
  `CATE` = rnd(point_estimate),
  `SE` = rnd(se_adj, 3),
  `95% CI low` = rnd(lower_95_adj), `95% CI high` = rnd(upper_95_adj)) %>%
  arrange(Slope, factor(Subgroup, levels = c("Above threshold","Below threshold","Difference (above - below)")))

# Table 4: methods comparison -- naive vs valid CI + variance components
t4 <- cmp %>% left_join(vcomp, by = "estimand") %>%
  transmute(
    Slope = slope_of(estimand), Subgroup = lab_map[estimand],
    `Naive CI width` = rnd(naive_width),
    `Boot-CV CI width` = rnd(bootcv_width),
    `Width ratio (valid/naive)` = rnd(width_ratio, 2),
    `Splitting SD (tau)` = rnd(sqrt(tau2_within)),
    `Sampling SE (sigma_BT)` = rnd(se_BT, 3),
    `Optimal B_CV` = rnd(optimal_B_CV, 0))

notes <- data.frame(Note = c(
  "Estimand: performance of the AUTOC UACR-threshold DISCOVERY procedure (trained on ~2788 patients), per Cai et al. 2025 -- not one fixed threshold value.",
  "Units: mL/min/1.73m2/yr. Positive CATE = intensive BP control preserves eGFR slope.",
  "Point estimates: 1000-iteration trial-stratified 60/40 Monte-Carlo CV (step 05).",
  "Confidence intervals: bootstrap-CV variance decomposition, Remark-3 adjusted (step 06b).",
  "Naive percentile interval (Table 4) is shown only to demonstrate it is NOT a valid CI: its width tracks splitting noise (tau), not sampling error (sigma_BT).",
  "chronic_above has optimal B_CV ~ 68 > the B_CV=10 used, so its SE is the least precise of the six; consider a B_CV=20 rerun before final submission."))

write_xlsx(list(
  `T1 full-data point`      = t1,
  `T2 CV reliability`       = t2,
  `T3 subgroup CATE (main)` = t3,
  `T4 methods comparison`   = t4,
  Notes                     = notes),
  "results/dissertation_tables.xlsx")

cat("Wrote:\n",
    " figures/fig1_forest_subgroup_cate.(pdf|png)\n",
    " figures/fig2_threshold_stability.(pdf|png)\n",
    " figures/fig3_variance_decomposition.(pdf|png)\n",
    " results/dissertation_tables.xlsx\n")
cat("\n--- T3 (headline) ---\n"); print(t3, row.names = FALSE)
