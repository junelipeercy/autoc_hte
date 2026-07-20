#!/bin/bash
#SBATCH --job-name=autoc_cv
#SBATCH --array=1-100%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=autoc_cv_%A_%a.log
# %A = array job id, %a = task number -> one log per task.
# 100 tasks x N_REP=10 = 1000 cross-validation iterations.

ml R

# --- ONE-TIME package install (uncomment for the very first run only) --------
# R -e 'install.packages(c("readxl","nlme","lspline","dplyr","writexl"), repos="https://cloud.r-project.org/")'

# N_REP / TRAIN_FRAC come in via --export (see below). Each task reads
# SLURM_ARRAY_TASK_ID / SLURM_CPUS_PER_TASK itself.
Rscript "05_autoc_cv.R"

# =============================================================================
# HOW TO SUBMIT  (run from ~, where analysis_long.xlsx + the .R files live)
#
# SMOKE TEST (2 tasks, 1 iteration each -> 2 CV iterations):
#   sbatch --array=1-2%2 --export=ALL,N_REP=1 run_autoc_cv_array.sbatch
#
# FULL RUN (100 tasks x 10 -> 1000 CV iterations):
#   sbatch --export=ALL,N_REP=10 run_autoc_cv_array.sbatch
#
# COMBINE after the array finishes (light; sh_dev node):
#   AUTOC_POINT=results/cv_point.csv Rscript 05b_combine_cv.R
#
# WATCH:   squeue -u $USER
# STATUS:  sacct -j <jobid> --format=JobID,State,ExitCode
# =============================================================================
