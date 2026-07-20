#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=autoc_boot_%x.log
# %x = job name; override it per cohort at submit time (see below).

ml R

# --- one-time package install (uncomment for the FIRST run only) -------------
# R -e 'install.packages(c("readxl","writexl","nlme","lspline","dplyr","tibble"), repos="https://cloud.r-project.org/")'

# Cohort is passed as the first argument after the script name.
# 03_autoc_bootstrap.R reads $SLURM_CPUS_PER_TASK to use all 8 cores.
Rscript 03_autoc_bootstrap.R "$1"

# -----------------------------------------------------------------------------
# Submit all five cohorts (each becomes its own concurrent job). --job-name
# sets %x, so each gets its own log (autoc_boot_POOLED.log, etc.):
#
#   sbatch --job-name=POOLED run_autoc_bootstrap.sbatch POOLED
#   sbatch --job-name=SPRINT run_autoc_bootstrap.sbatch SPRINT
#   sbatch --job-name=ACCORD run_autoc_bootstrap.sbatch ACCORD
#   sbatch --job-name=AASK   run_autoc_bootstrap.sbatch AASK
#   sbatch --job-name=MDRD   run_autoc_bootstrap.sbatch MDRD
#
# Watch:   squeue -u $USER
# Stream:  tail -f autoc_boot_POOLED.log
# -----------------------------------------------------------------------------
