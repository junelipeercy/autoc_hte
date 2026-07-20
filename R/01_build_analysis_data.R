# =============================================================================
# 01_build_analysis_data.R  --  Build ONE master analysis dataset.
#
# Combines the four cohorts (SPRINT, ACCORD, AASK, MDRD) into a single long
# dataset that already carries all three baseline covariates (egfr0, uacr,
# kfrs), so downstream scripts never have to piece datasets together again.
#
# Inputs  (data/raw/, git-ignored, IRB-private):
#   combined_long.xlsx   SPRINT + ACCORD (trial = first char of id: S / A)
#   aask_dataset.xlsx    AASK
#   mdrd_dataset.xlsx    MDRD
#
# Output:
#   data/derived/analysis_long.xlsx
#     one row per eGFR visit; columns: id, trial, arm, egfr, days, egfr0,
#     uacr, kfrs.  egfr0 / uacr / kfrs are constant within id (true baseline).
#
# Assumption: arm is coded consistently across cohorts (1 = intensive BP arm).
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
})

raw <- "data/raw"
keep_cols <- c("id", "trial", "arm", "egfr", "days", "egfr0", "uacr", "kfrs")

# --- SPRINT + ACCORD ---------------------------------------------------------
sa <- read_excel(file.path(raw, "combined_long.xlsx"))
names(sa) <- tolower(names(sa))
sa <- sa %>%
  mutate(trial = ifelse(substr(id, 1, 1) == "S", "SPRINT", "ACCORD"),
         id    = as.character(id)) %>%
  select(any_of(keep_cols))

# --- AASK  (prefix ids with "AA" to keep them unique when pooled) ------------
aask <- read_excel(file.path(raw, "aask_dataset.xlsx"))
names(aask) <- tolower(names(aask))
aask <- aask %>%
  mutate(trial = "AASK",
         id    = paste0("AA", id)) %>%
  select(any_of(keep_cols))

# --- MDRD  (prefix ids with "M") ---------------------------------------------
mdrd <- read_excel(file.path(raw, "mdrd_dataset.xlsx"))
names(mdrd) <- tolower(names(mdrd))
mdrd <- mdrd %>%
  mutate(trial = "MDRD",
         id    = paste0("M", id)) %>%
  select(any_of(keep_cols))

# --- Pool + clean ------------------------------------------------------------
analysis_long <- bind_rows(sa, aask, mdrd) %>%
  filter(!is.na(egfr), !is.na(days), !is.na(arm), !is.na(id))

# --- Report so you can eyeball it before saving ------------------------------
cat("\nRows per trial:\n")
print(count(analysis_long, trial))

cat("\nUnique patients per trial:\n")
print(analysis_long %>% distinct(id, trial) %>% count(trial))

cat("\nBaseline covariate completeness (patients with non-missing value):\n")
print(
  analysis_long %>%
    distinct(id, trial, egfr0, uacr, kfrs) %>%
    group_by(trial) %>%
    summarise(n          = n(),
              egfr0_ok    = sum(!is.na(egfr0)),
              uacr_ok     = sum(!is.na(uacr)),
              kfrs_ok     = sum(!is.na(kfrs)),
              .groups = "drop")
)

# --- Save --------------------------------------------------------------------
dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
write_xlsx(analysis_long, "data/derived/analysis_long.xlsx")
cat("\nWrote data/derived/analysis_long.xlsx  (",
    nrow(analysis_long), "rows )\n")
