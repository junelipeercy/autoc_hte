# =============================================================================
# 03_autoc_bootstrap.R  --  ONE array task's share of the chronic+total AUTOC
# bootstrap, for ANY cohort (POOLED or a single trial). Self-contained.
#
# Same validated logic as the pooled run; the only additions are a COHORT switch
# (subset the data) and COHORT-tagged fragment filenames so trials don't collide.
#
# Run as part of a Slurm job array (run_autoc_bootstrap_array.sbatch). Each task:
#   - reads its number from  SLURM_ARRAY_TASK_ID  (1..100)
#   - runs N_BOOT bootstraps (default 10) with seed = SEED_BASE + task_id
#   - runs them in parallel across SLURM_CPUS_PER_TASK cores
#   - writes two COHORT-tagged fragment files to the output dir:
#         frag_scalars_<COHORT>_<task>.csv   one row per bootstrap
#         frag_curves_<COHORT>_<task>.csv    long AUTOC curves (3 bm x n_pct)
#
# Env vars (all optional):
#   COHORT                POOLED | SPRINT | ACCORD | AASK | MDRD   (default POOLED)
#   SLURM_ARRAY_TASK_ID   which task (default 1)
#   SLURM_CPUS_PER_TASK   cores for mclapply (default 4)
#   N_BOOT                bootstraps this task runs (default 10; set 1 to smoke-test)
#   AUTOC_DATA            path to analysis_long.xlsx (default ~/analysis_long.xlsx)
#   AUTOC_OUT             output dir for fragments (default $SCRATCH, else ".")
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(nlme); library(lspline)
  library(dplyr); library(parallel)
})

# ---- config -----------------------------------------------------------------
COHORT  <- toupper(Sys.getenv("COHORT", "POOLED"))
TASK    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
N_BOOT  <- as.integer(Sys.getenv("N_BOOT", "10"))
SEED_BASE <- 20260714L
KNOT    <- 4
TARGET_CHRONIC <- 1.0     # UACR threshold: chronic CATE crosses +1.0
TARGET_TOTAL   <- 0.75    # UACR threshold: total   CATE crosses +0.75

# Percentile depth: smaller trials capped at 85 (less noisy extremes).
NPCT_MAP <- c(POOLED = 90, SPRINT = 90, ACCORD = 85, AASK = 90, MDRD = 85)
if (!COHORT %in% names(NPCT_MAP)) stop("Unknown COHORT: ", COHORT)
N_PCT <- NPCT_MAP[[COHORT]]

DATA <- Sys.getenv("AUTOC_DATA", unset = file.path(Sys.getenv("HOME"), "analysis_long.xlsx"))
OUT  <- Sys.getenv("AUTOC_OUT",  unset = Sys.getenv("SCRATCH", unset = "."))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("COHORT=%s | Task %d | N_BOOT=%d | cores=%d | knot=%d | n_pct=%d\nData: %s\nOut:  %s\n",
            COHORT, TASK, N_BOOT, NCORES, KNOT, N_PCT, DATA, OUT))

# ---- functions --------------------------------------------------------------
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

# AUTOC curve for one biomarker: CATE (both slopes) at each percentile cutoff
autoc_curve <- function(df, bm, knot = KNOT, n_pct = N_PCT) {
  base <- df[!duplicated(df$id), ]
  v <- switch(bm, egfr = base$egfr0, uacr = base$uacr, kfre = base$kfrs)
  v <- sort(unique(v[!is.na(v)]))
  if (length(v) < n_pct) return(NULL)
  pct <- quantile(v, probs = seq(0.01, 1, by = 0.01), na.rm = TRUE)
  cutoffs <- if (bm == "egfr") sort(pct, decreasing = TRUE)[1:n_pct] else sort(pct)[1:n_pct]
  tot <- chr <- rep(NA_real_, n_pct)
  for (j in seq_along(cutoffs)) {
    i <- cutoffs[j]
    sub <- switch(bm, egfr = df[df$egfr0 <= i, ], uacr = df[df$uacr >= i, ], kfre = df[df$kfrs >= i, ])
    sub <- sub[!is.na(sub$egfr) & !is.na(sub$days) & !is.na(sub$arm), ]
    a <- tryCatch(fit_ate(sub, knot), error = function(e) c(total = NA_real_, chronic = NA_real_))
    tot[j] <- a["total"]; chr[j] <- a["chronic"]
  }
  data.frame(biomarker = bm, percentile = seq_len(n_pct),
             cov_value = as.numeric(cutoffs), total_cate = tot, chronic_cate = chr)
}

cross_threshold <- function(curve, col, target) {
  y <- curve[[col]]; x <- curve$cov_value
  idx <- which(y >= target)
  if (length(idx) == 0) return(c(thr = NA_real_, crossed = 0))
  j <- idx[1]
  if (j == 1) return(c(thr = x[1], crossed = 1))
  a0 <- y[j - 1]; a1 <- y[j]; c0 <- x[j - 1]; c1 <- x[j]
  thr <- if (is.na(a0) || is.na(a1) || a1 == a0) x[j] else c0 + (target - a0) / (a1 - a0) * (c1 - c0)
  c(thr = thr, crossed = 1)
}

subgroup_ate <- function(df, thr, slope) {
  ab <- df[df$uacr >= thr, ]; ab <- ab[!is.na(ab$egfr) & !is.na(ab$days) & !is.na(ab$arm), ]
  be <- df[df$uacr <  thr, ]; be <- be[!is.na(be$egfr) & !is.na(be$days) & !is.na(be$arm), ]
  aa <- tryCatch(fit_ate(ab), error = function(e) c(total = NA, chronic = NA))
  bb <- tryCatch(fit_ate(be), error = function(e) c(total = NA, chronic = NA))
  c(above = unname(aa[slope]), below = unname(bb[slope]),
    diff = unname(aa[slope] - bb[slope]),
    n_above = length(unique(ab$id)), n_below = length(unique(be$id)))
}

# ---- load, subset to COHORT, pre-split by patient ---------------------------
al <- read_excel(DATA)
al <- al[!is.na(al$egfr) & !is.na(al$days) & !is.na(al$arm) & !is.na(al$id), ]
if (COHORT != "POOLED") al <- al[al$trial == COHORT, ]
if (nrow(al) == 0) stop("No rows for COHORT=", COHORT, " (check the 'trial' column).")
id_rows    <- split(al, al$id)
unique_ids <- names(id_rows)
cat(COHORT, "patients:", length(unique_ids), " rows:", nrow(al), "\n")

# ---- one bootstrap ----------------------------------------------------------
one_boot <- function(b) {
  gid <- (TASK - 1L) * N_BOOT + b
  sampled <- sample(unique_ids, length(unique_ids), replace = TRUE)
  df_boot <- bind_rows(Map(function(sid, i) { t <- id_rows[[sid]]; t$id <- paste0("b", i); t },
                           sampled, seq_along(sampled)))

  ce <- autoc_curve(df_boot, "egfr")
  cu <- autoc_curve(df_boot, "uacr")
  ck <- autoc_curve(df_boot, "kfre")
  if (is.null(ce) || is.null(cu) || is.null(ck)) return(NULL)

  A <- function(cate) auc_trapezoid(N_PCT:1, cate - cate[1])
  a_e_c <- A(ce$chronic_cate); a_u_c <- A(cu$chronic_cate); a_k_c <- A(ck$chronic_cate)
  a_e_t <- A(ce$total_cate);   a_u_t <- A(cu$total_cate);   a_k_t <- A(ck$total_cate)

  overall_chronic <- cu$chronic_cate[1]
  overall_total   <- cu$total_cate[1]

  labs <- c("eGFR", "UACR", "KFRE")
  winner_chronic <- labs[which.max(c(a_e_c, a_u_c, a_k_c))]
  winner_total   <- labs[which.max(c(a_e_t, a_u_t, a_k_t))]

  tc <- cross_threshold(cu, "chronic_cate", TARGET_CHRONIC)
  tt <- cross_threshold(cu, "total_cate",   TARGET_TOTAL)
  sc <- if (tc["crossed"] == 1) subgroup_ate(df_boot, tc["thr"], "chronic")
        else c(above = NA, below = NA, diff = NA, n_above = NA, n_below = NA)
  st <- if (tt["crossed"] == 1) subgroup_ate(df_boot, tt["thr"], "total")
        else c(above = NA, below = NA, diff = NA, n_above = NA, n_below = NA)

  n_na <- sum(is.na(c(ce$chronic_cate, cu$chronic_cate, ck$chronic_cate,
                      ce$total_cate,  cu$total_cate,  ck$total_cate)))

  scalars <- data.frame(
    cohort = COHORT, task = TASK, boot = gid, seed = SEED_BASE + TASK,
    autoc_chronic_egfr = a_e_c, autoc_chronic_uacr = a_u_c, autoc_chronic_kfre = a_k_c,
    autoc_total_egfr = a_e_t, autoc_total_uacr = a_u_t, autoc_total_kfre = a_k_t,
    overall_chronic = overall_chronic, overall_total = overall_total,
    winner_chronic = winner_chronic, winner_total = winner_total,
    thr_chronic = unname(tc["thr"]), crossed_chronic = unname(tc["crossed"]),
    thr_total   = unname(tt["thr"]), crossed_total   = unname(tt["crossed"]),
    sg_chronic_above = sc["above"], sg_chronic_below = sc["below"], sg_chronic_diff = sc["diff"],
    n_chronic_above = sc["n_above"], n_chronic_below = sc["n_below"],
    sg_total_above = st["above"], sg_total_below = st["below"], sg_total_diff = st["diff"],
    n_total_above = st["n_above"], n_total_below = st["n_below"],
    n_na_fits = n_na, stringsAsFactors = FALSE, row.names = NULL)

  curves <- rbind(ce, cu, ck); curves$boot <- gid; curves$cohort <- COHORT
  list(scalars = scalars, curves = curves)
}

# ---- run this task's bootstraps in parallel ---------------------------------
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED_BASE + TASK)
t0 <- Sys.time()
res <- mclapply(seq_len(N_BOOT), one_boot, mc.cores = NCORES, mc.set.seed = TRUE)
res <- Filter(Negate(is.null), res)
cat(COHORT, "task", TASK, "done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

scalars <- do.call(rbind, lapply(res, `[[`, "scalars"))
curves  <- do.call(rbind, lapply(res, `[[`, "curves"))

tag <- sprintf("%s_%03d", COHORT, TASK)
write.csv(scalars, file.path(OUT, paste0("frag_scalars_", tag, ".csv")), row.names = FALSE)
write.csv(curves,  file.path(OUT, paste0("frag_curves_",  tag, ".csv")), row.names = FALSE)
cat("Wrote fragments to", OUT, "\n")
