# =============================================================================
# 06_autoc_bootcv.R  --  Bootstrap-cross-validation (Cai, ..., Tian 2025,
# "Bootstrapping the Cross-Validation Estimate", Ann. Appl. Stat.; arXiv:2307.00260)
# applied to the UACR-threshold AUTOC pipeline.
#
# WHY: repeated-split Monte-Carlo CV (05_autoc_cv.R) gives a good POINT estimate
# but its spread across splits measures only the train/test SPLITTING noise
# (tau^2), NOT sampling uncertainty. Percentile CIs from it are therefore not
# valid CIs and get NARROWER the more splits you run. This script adds the outer
# bootstrap loop that supplies the missing sampling variance (sigma^2_BT).
#
# THEIR ALGORITHM 2, as implemented here:
#   for b in 1..B_BOOT:
#     - draw patient-level bootstrap weights W_b ~ Multinomial(n, 1/n)
#     - for k in 1..B_CV:
#         * split the ORIGINAL patients into train (m_adj) / test (n - m_adj),
#           trial-stratified                       <- SPLIT FIRST ...
#         * apply W_b to each side separately      <- ... THEN WEIGHT
#           (train and test stay disjoint => no leakage, unlike naive bootstrap
#           which resamples first and lets duplicates straddle both sides)
#         * run the full pipeline: AUTOC threshold search on the weighted train,
#           subgroup CATEs on the weighted test  ->  theta*_bk
#
# m_adj: a bootstrapped training set holds only 0.632*m DISTINCT patients, so
# the training size is inflated by minimising their loss (eq. 2) with
# lambda0 = 0.368, trading effective training size against test-set shrinkage.
#
# The variance decomposition theta*_bk = theta0 + eps*_b + eps_bk is fitted in
# 06b_combine_bootcv.R (closed-form moment estimators from their step 11).
#
# Weighting note: for an LMM with per-patient random effects, "weight patient i
# by W_i" is operationalised as REPLICATING that patient W_i times with distinct
# ids, so each replicate gets its own random effect -- same convention as
# 03_autoc_bootstrap.R.
#
# Remark 4 of the paper covers our data structure: longitudinal measurements are
# correlated within patient, but patients are i.i.d., so both the CV split and
# the bootstrap are done at the PATIENT level.
#
# Env vars (all optional):
#   SLURM_ARRAY_TASK_ID   which task (default 1)
#   SLURM_CPUS_PER_TASK   cores for mclapply (default 4)
#   N_BOOT_TASK           bootstrap draws THIS task runs (default 4)
#   B_CV                  inner CV splits per bootstrap (default 10; paper 10-20)
#   TRAIN_FRAC            nominal training fraction m/n (default 0.60)
#   AUTOC_DATA            path to analysis_long.xlsx (default ~/analysis_long.xlsx)
#   AUTOC_OUT             output dir for fragments (default $SCRATCH, else ".")
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(nlme); library(lspline); library(dplyr); library(parallel)
})

# ---- config -----------------------------------------------------------------
TASK        <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
NCORES      <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
N_BOOT_TASK <- as.integer(Sys.getenv("N_BOOT_TASK", "4"))
B_CV        <- as.integer(Sys.getenv("B_CV", "10"))
TRAIN_FRAC  <- as.numeric(Sys.getenv("TRAIN_FRAC", "0.60"))
SEED_BASE   <- 20260721L
KNOT        <- 4
N_PCT       <- 90
TARGET_CHRONIC <- 1.00
TARGET_TOTAL   <- 0.75
LAMBDA0     <- 0.368        # = 1 - 0.632, the paper's recommended penalty

DATA <- Sys.getenv("AUTOC_DATA", unset = file.path(Sys.getenv("HOME"), "analysis_long.xlsx"))
OUT  <- Sys.getenv("AUTOC_OUT",  unset = Sys.getenv("SCRATCH", unset = "."))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- engine (identical math to 03/05) ---------------------------------------
auc_trapezoid <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

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

uacr_curve <- function(df, knot = KNOT, n_pct = N_PCT) {
  base <- df[!duplicated(df$id), ]
  v <- sort(unique(base$uacr[!is.na(base$uacr)]))
  if (length(v) < n_pct) return(NULL)
  pct     <- quantile(v, probs = seq(0.01, 1, by = 0.01), na.rm = TRUE)
  cutoffs <- sort(pct)[1:n_pct]
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

cross_threshold <- function(curve, col, target) {
  y <- curve[[col]]; x <- curve$cov_value
  idx <- which(y >= target)
  if (length(idx) == 0) return(c(thr = NA, crossed = 0))
  j <- idx[1]
  if (j == 1) return(c(thr = x[1], crossed = 1))
  a0 <- y[j - 1]; a1 <- y[j]; c0 <- x[j - 1]; c1 <- x[j]
  thr <- if (is.na(a0) || is.na(a1) || a1 == a0) x[j]
         else c0 + (target - a0) / (a1 - a0) * (c1 - c0)
  c(thr = thr, crossed = 1)
}

NA_SG <- c(above = NA, below = NA, diff = NA, n_above = NA, n_below = NA)
eval_subgroups <- function(df, thr, slope) {
  ab <- df[df$uacr >= thr, ]; ab <- ab[!is.na(ab$egfr) & !is.na(ab$days) & !is.na(ab$arm), ]
  be <- df[df$uacr <  thr, ]; be <- be[!is.na(be$egfr) & !is.na(be$days) & !is.na(be$arm), ]
  aa <- tryCatch(fit_ate(ab), error = function(e) c(total = NA, chronic = NA))
  bb <- tryCatch(fit_ate(be), error = function(e) c(total = NA, chronic = NA))
  c(above = unname(aa[slope]), below = unname(bb[slope]),
    diff = unname(aa[slope] - bb[slope]),
    n_above = length(unique(ab$id)), n_below = length(unique(be$id)))
}

# ---- load, index patients ---------------------------------------------------
al <- read_excel(DATA)
al <- al[!is.na(al$egfr) & !is.na(al$days) & !is.na(al$arm) & !is.na(al$id) & !is.na(al$uacr), ]
if (nrow(al) == 0) stop("No usable rows in ", DATA)
id_rows      <- split(al, al$id)
id_trial     <- al[!duplicated(al$id), c("id", "trial")]
ids_by_trial <- split(id_trial$id, id_trial$trial)
ALL_IDS      <- id_trial$id
N_PT         <- length(ALL_IDS)

# ---- m_adj: paper eq. (2) ---------------------------------------------------
# minimise (m_adj/(m/0.632) - 1)^2 + lambda0*((n-m)/(n-m_adj) - 1)^2
M_NOM <- floor(TRAIN_FRAC * N_PT)
madj_loss <- function(a) {
  (a / (M_NOM / 0.632) - 1)^2 + LAMBDA0 * ((N_PT - M_NOM) / (N_PT - a) - 1)^2
}
M_ADJ <- round(optimize(madj_loss, interval = c(M_NOM, N_PT - 50))$minimum)
ADJ_FRAC <- M_ADJ / N_PT

cat(sprintf(paste0("BOOT-CV | Task %d | bootstraps=%d | B_CV=%d | cores=%d\n",
                   "n=%d patients | m(nominal)=%d (%.0f%%) | m_adj=%d (%.1f%%) | test=%d\n",
                   "Data: %s\nOut:  %s\n"),
            TASK, N_BOOT_TASK, B_CV, NCORES, N_PT, M_NOM, 100 * TRAIN_FRAC,
            M_ADJ, 100 * ADJ_FRAC, N_PT - M_ADJ, DATA, OUT))
cat("By trial:", paste(names(ids_by_trial), lengths(ids_by_trial), sep = "=", collapse = "  "), "\n")

# ---- helpers ----------------------------------------------------------------
# trial-stratified split of the ORIGINAL patients at the m_adj fraction
stratified_split <- function() {
  train <- unlist(lapply(ids_by_trial, function(ids)
    sample(ids, floor(ADJ_FRAC * length(ids)))), use.names = FALSE)
  list(train = train, test = setdiff(ALL_IDS, train))
}

# apply bootstrap weights to one side of the split: replicate patient i W_i
# times with distinct ids (W_i = 0 patients simply drop out).
expand_by_weights <- function(ids, W) {
  sel <- rep(ids, times = as.integer(W[ids]))
  if (length(sel) == 0) return(NULL)
  bind_rows(Map(function(sid, j) { t <- id_rows[[sid]]; t$id <- paste0("w", j); t },
                sel, seq_along(sel)))
}

# ---- one (b, k) unit: one bootstrap weight vector x one CV split -------------
one_unit <- function(b_local, k, W, b_global) {
  sp <- stratified_split()                       # SPLIT FIRST (original data)
  tr <- expand_by_weights(sp$train, W)           # ... THEN WEIGHT, disjointly
  te <- expand_by_weights(sp$test,  W)
  if (is.null(tr) || is.null(te)) return(NULL)

  cu <- uacr_curve(tr)
  if (is.null(cu)) return(NULL)

  tc <- cross_threshold(cu, "chronic_cate", TARGET_CHRONIC)
  tt <- cross_threshold(cu, "total_cate",   TARGET_TOTAL)
  sc <- if (tc["crossed"] == 1) eval_subgroups(te, tc["thr"], "chronic") else NA_SG
  st <- if (tt["crossed"] == 1) eval_subgroups(te, tt["thr"], "total")   else NA_SG

  data.frame(
    boot = b_global, cv = k, task = TASK,
    m_adj = M_ADJ, n_train_w = length(unique(tr$id)), n_test_w = length(unique(te$id)),
    autoc_chronic_train = auc_trapezoid(N_PCT:1, cu$chronic_cate - cu$chronic_cate[1]),
    autoc_total_train   = auc_trapezoid(N_PCT:1, cu$total_cate   - cu$total_cate[1]),
    thr_chronic = unname(tc["thr"]), crossed_chronic = unname(tc["crossed"]),
    thr_total   = unname(tt["thr"]), crossed_total   = unname(tt["crossed"]),
    # theta*_bk : the four estimands of interest
    chronic_above = sc["above"], chronic_below = sc["below"], chronic_diff = sc["diff"],
    total_above   = st["above"], total_below   = st["below"], total_diff   = st["diff"],
    n_chronic_above = sc["n_above"], n_chronic_below = sc["n_below"],
    n_total_above   = st["n_above"], n_total_below   = st["n_below"],
    stringsAsFactors = FALSE, row.names = NULL)
}

# ---- run: pre-draw this task's bootstrap weights, then parallelise (b,k) ----
# Weights are drawn in the PARENT so all B_CV splits within a bootstrap share
# the same W_b (required by the variance decomposition).
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED_BASE + TASK)

boot_ids <- (TASK - 1L) * N_BOOT_TASK + seq_len(N_BOOT_TASK)
W_list <- lapply(seq_len(N_BOOT_TASK), function(i) {
  w <- as.vector(rmultinom(1, size = N_PT, prob = rep(1 / N_PT, N_PT)))
  setNames(w, ALL_IDS)
})

grid <- expand.grid(b_local = seq_len(N_BOOT_TASK), k = seq_len(B_CV))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(i) {
  bl <- grid$b_local[i]; k <- grid$k[i]
  tryCatch(one_unit(bl, k, W_list[[bl]], boot_ids[bl]), error = function(e) NULL)
}, mc.cores = NCORES, mc.set.seed = TRUE)
res <- Filter(Negate(is.null), res)
cat("Task", TASK, "done:", length(res), "of", nrow(grid), "units in",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

scalars <- do.call(rbind, res)
tag <- sprintf("%03d", TASK)
write.csv(scalars, file.path(OUT, paste0("frag_bootcv_", tag, ".csv")), row.names = FALSE)
cat("Wrote frag_bootcv_", tag, ".csv to ", OUT, "\n", sep = "")
