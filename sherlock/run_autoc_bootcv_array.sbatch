#!/bin/bash
#SBATCH --job-name=autoc_bootcv
#SBATCH --array=1-100%25
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=autoc_bootcv_%A_%a.log
# %A = array job id, %a = task number -> one log per task.
#
# 100 tasks x N_BOOT_TASK=4 bootstraps x B_CV=10 splits = 4000 pipeline runs
#   (= B_BOOT 400, the paper's recommended setting).
# Each task runs 40 units across 8 cores.

ml R

# --- ONE-TIME package install (uncomment for the very first run only) --------
# R -e 'install.packages(c("readxl","nlme","lspline","dplyr","writexl"), repos="https://cloud.r-project.org/")'

Rscript "06_autoc_bootcv.R"

# =============================================================================
# HOW TO SUBMIT  (run from ~, where analysis_long.xlsx + the .R files live)
#
# SMOKE TEST (1 task, 1 bootstrap x 2 splits -> 2 units; confirm it runs):
#   sbatch --array=1-1 --export=ALL,N_BOOT_TASK=1,B_CV=2 run_autoc_bootcv_array.sbatch
#
# FULL RUN (B_BOOT=400, B_CV=10):
#   sbatch --export=ALL,N_BOOT_TASK=4,B_CV=10 run_autoc_bootcv_array.sbatch
#
# If the combine step reports optimal_B_CV >> 10, rerun with a larger B_CV:
#   sbatch --export=ALL,N_BOOT_TASK=4,B_CV=20 run_autoc_bootcv_array.sbatch
#
# COMBINE after the array finishes (light; sh_dev node). Supply the MC-CV raw
# file so the point estimate is the correct Err^CV_m and you get the
# naive-vs-valid CI comparison table:
#   N_PATIENTS=4647 AUTOC_MCCV=~/autoc_cv_raw.csv Rscript 06b_combine_bootcv.R
#
# WATCH:   squeue -u $USER
# STATUS:  sacct -j <jobid> --format=JobID,State,Elapsed,ExitCode
# =============================================================================
