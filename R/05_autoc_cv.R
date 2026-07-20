# =============================================================================
# 05_autoc_cv.R  --  PART 2: cross-validating the UACR-threshold DISCOVERY.
#
# One array task's share of a 1000-iteration Monte-Carlo cross-validation on the
# POOLED cohort. UACR is fixed as the ranking covariate (already chosen as the
# best HTE marker by the step-3 bootstrap); here we only ask whether the AUTOC
# threshold-finding PROCESS gives a reliable, roughly unbiased cutoff.
#
# Each CV iteration:
#   1. SPLIT patients 60/40 (TRAIN/TEST), WITHOUT replacement, STRATIFIED by
#      trial so every split keeps the pooled SPRINT/ACCORD/AASK/MDRD mix. The
#      split is by patient id -> a patient's whole trajectory stays on one side
#      (no leakage).
#   2. TRAIN (60%): build the UACR-ranked CATE curve and find the threshold
#      where chronic CATE crosses +1.00 and where total CATE crosses +0.75.
#      Also record the in-sample above/below effect at that threshold (this is
#      the OPTIMISTIC number; comparing it to TEST measures the method's bias).
#   3. TEST (40%): split the held-out patients at the TRAIN-derived threshold
#      and estimate the treatment effect above vs below -- the HONEST estimate.
#   4. If a slope's target is never crossed in TRAIN, that slope's iteration is
#      marked crossed=0 (kept for the failure-rate denominator, dropped from the
#      threshold/effect summaries downstream).
#
# Mirrors 03_autoc_bootstrap.R's engine (fit_ate / cross_threshold / subgroup
# fits); the only new machinery is the stratified split and TRAIN-vs-TEST split.
#
# Env vars (all optional):
#   SLURM_ARRAY_TASK_ID   which task (default 1)
#   SLURM_CPUS_PER_TASK   cores for mclapply (default 4)
#   N_REP                 CV iterations this task runs (default 10; 1 to smoke-test)
#   TRAIN_FRAC            training fraction (default 0.60)
#   AUTOC_DATA            path to analysis_long.xlsx (default ~/analysis_long.xlsx)
#   AUTOC_OUT             output dir for fragments (default $SCRATCH, else ".")
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(nlme); library(lspline); library(dplyr); library(parallel)
})

# ---- config -----------------------------------------------------------------
TASK       <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
NCORES     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
N_REP      <- as.integer(Sys.getenv("N_REP", "10"))
TRAIN_FRAC <- as.numeric(Sys.getenv("TRAIN_FRAC", "0.60"))
SEED_BASE  <- 20260720L
KNOT       <- 4
N_PCT      <- 90            # pooled percentile depth (matches step 3 POOLED)
TARGET_CHRONIC <- 1.00      # UACR threshold: chronic CATE crosses +1.00
TARGET_TOTAL   <- 0.75      # UACR threshold: total   CATE crosses +0.75

DATA <- Sys.getenv("AUTOC_DATA", unset = file.path(Sys.getenv("HOME"), "analysis_long.xlsx"))
OUT  <- Sys.getenv("AUTOC_OUT",  unset = Sys.getenv("SCRATCH", unset = "."))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("CV | Task %d | N_REP=%d | train_frac=%.2f | cores=%d | knot=%d | n_pct=%d\nData: %s\nOut:  %s\n",
            TASK, N_REP, TRAIN_FRAC, NCORES, KNOT, N_PCT, DATA, OUT))

# ---- engine (identical math to 03_autoc_bootstrap.R) ------------------------
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

# UACR-ranked CATE curve on a dataset (TRAIN fold), both slopes.
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

# above/below treatment effect for one slope, + subgroup sizes + convergence.
NA_SG <- c(above = NA, below = NA, diff = NA, n_above = NA, n_below = NA,
           ok_above = 0, ok_below = 0)
eval_subgroups <- function(df, thr, slope) {
  ab <- df[df$uacr >= thr, ]; ab <- ab[!is.na(ab$egfr) & !is.na(ab$days) & !is.na(ab$arm), ]
  be <- df[df$uacr <  thr, ]; be <- be[!is.na(be$egfr) & !is.na(be$days) & !is.na(be$arm), ]
  aa <- tryCatch(fit_ate(ab), error = function(e) c(total = NA, chronic = NA))
  bb <- tryCatch(fit_ate(be), error = function(e) c(total = NA, chronic = NA))
  c(above = unname(aa[slope]), below = unname(bb[slope]),
    diff = unname(aa[slope] - bb[slope]),
    n_above = length(unique(ab$id)), n_below = length(unique(be$id)),
    ok_above = as.numeric(!is.na(aa[slope])), ok_below = as.numeric(!is.na(bb[slope])))
}

# ---- load pooled data, index patients by trial for stratified splitting -----
al <- read_excel(DATA)
al <- al[!is.na(al$egfr) & !is.na(al$days) & !is.na(al$arm) & !is.na(al$id) & !is.na(al$uacr), ]
if (nrow(al) == 0) stop("No usable rows in ", DATA)
id_trial   <- al[!duplicated(al$id), c("id", "trial")]
ids_by_trial <- split(id_trial$id, id_trial$trial)
cat("Pooled patients:", nrow(id_trial), " rows:", nrow(al), "\n")
cat("By trial:", paste(names(ids_by_trial), lengths(ids_by_trial), sep = "=", collapse = "  "), "\n")

# trial-stratified 60/40 split of patient ids (no replacement).
stratified_split <- function() {
  train <- unlist(lapply(ids_by_trial, function(ids) {
    n <- length(ids); sample(ids, floor(TRAIN_FRAC * n))
  }), use.names = FALSE)
  list(train = train, test = setdiff(id_trial$id, train))
}

# ---- one CV iteration -------------------------------------------------------
one_rep <- function(r) {
  gid <- (TASK - 1L) * N_REP + r
  sp  <- stratified_split()
  tr  <- al[al$id %in% sp$train, ]
  te  <- al[al$id %in% sp$test,  ]

  cu <- uacr_curve(tr)                       # TRAIN UACR curve
  if (is.null(cu)) return(NULL)              # too few distinct UACR (won't happen pooled)

  autoc_chronic <- auc_trapezoid(N_PCT:1, cu$chronic_cate - cu$chronic_cate[1])
  autoc_total   <- auc_trapezoid(N_PCT:1, cu$total_cate   - cu$total_cate[1])
  overall_chronic_train <- cu$chronic_cate[1]
  overall_total_train   <- cu$total_cate[1]

  tc <- cross_threshold(cu, "chronic_cate", TARGET_CHRONIC)   # threshold search
  tt <- cross_threshold(cu, "total_cate",   TARGET_TOTAL)

  # in-sample (TRAIN) effect at that threshold  -> optimism reference
  sc_tr <- if (tc["crossed"] == 1) eval_subgroups(tr, tc["thr"], "chronic") else NA_SG
  st_tr <- if (tt["crossed"] == 1) eval_subgroups(tr, tt["thr"], "total")   else NA_SG
  # held-out (TEST) effect at that threshold   -> honest estimate
  sc_te <- if (tc["crossed"] == 1) eval_subgroups(te, tc["thr"], "chronic") else NA_SG
  st_te <- if (tt["crossed"] == 1) eval_subgroups(te, tt["thr"], "total")   else NA_SG

  # whole-fold overall ATE (context for the subgroup numbers)
  ov_te <- tryCatch(fit_ate(te), error = function(e) c(total = NA, chronic = NA))

  data.frame(
    rep = gid, task = TASK, seed = SEED_BASE + TASK,
    n_train = length(sp$train), n_test = length(sp$test),
    autoc_chronic_train = autoc_chronic, autoc_total_train = autoc_total,
    overall_chronic_train = overall_chronic_train, overall_total_train = overall_total_train,
    overall_chronic_test = unname(ov_te["chronic"]), overall_total_test = unname(ov_te["total"]),

    # ---- chronic slope (target CATE > 1.00) ----
    thr_chronic = unname(tc["thr"]), crossed_chronic = unname(tc["crossed"]),
    boundary_chronic = unname(tc["boundary"]), ncross_chronic = unname(tc["ncross"]),
    tr_chronic_above = sc_tr["above"], tr_chronic_below = sc_tr["below"], tr_chronic_diff = sc_tr["diff"],
    te_chronic_above = sc_te["above"], te_chronic_below = sc_te["below"], te_chronic_diff = sc_te["diff"],
    n_te_chronic_above = sc_te["n_above"], n_te_chronic_below = sc_te["n_below"],
    ok_te_chronic_above = sc_te["ok_above"], ok_te_chronic_below = sc_te["ok_below"],

    # ---- total slope (target CATE > 0.75) ----
    thr_total = unname(tt["thr"]), crossed_total = unname(tt["crossed"]),
    boundary_total = unname(tt["boundary"]), ncross_total = unname(tt["ncross"]),
    tr_total_above = st_tr["above"], tr_total_below = st_tr["below"], tr_total_diff = st_tr["diff"],
    te_total_above = st_te["above"], te_total_below = st_te["below"], te_total_diff = st_te["diff"],
    n_te_total_above = st_te["n_above"], n_te_total_below = st_te["n_below"],
    ok_te_total_above = st_te["ok_above"], ok_te_total_below = st_te["ok_below"],

    stringsAsFactors = FALSE, row.names = NULL)
}

# ---- run this task's CV iterations in parallel ------------------------------
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED_BASE + TASK)
t0 <- Sys.time()
res <- mclapply(seq_len(N_REP), one_rep, mc.cores = NCORES, mc.set.seed = TRUE)
res <- Filter(Negate(is.null), res)
cat("CV task", TASK, "done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

scalars <- do.call(rbind, res)
tag <- sprintf("%03d", TASK)
write.csv(scalars, file.path(OUT, paste0("frag_cv_", tag, ".csv")), row.names = FALSE)
cat("Wrote fragment frag_cv_", tag, ".csv to ", OUT, "\n", sep = "")
