#!/bin/bash
#SBATCH --job-name=autoc_boot
#SBATCH --array=1-100%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=autoc_boot_%A_%a.log
# %A = array job id, %a = task number -> one log per task.

ml R

# --- ONE-TIME package install (uncomment for the very first run only) --------
# R -e 'install.packages(c("readxl","nlme","lspline","dplyr","tidyr","writexl"), repos="https://cloud.r-project.org/")'

# COHORT and N_BOOT come in via --export (see below). Each task reads
# SLURM_ARRAY_TASK_ID / SLURM_CPUS_PER_TASK itself.
Rscript "03_autoc_bootstrap.R"

# =============================================================================
# HOW TO SUBMIT  (run from ~, where analysis_long.xlsx + the .R files live)
#
# SMOKE TEST one cohort (2 tasks, 1 boot each):
#   sbatch --array=1-2%2 --export=ALL,N_BOOT=1,COHORT=ACCORD run_autoc_bootstrap_array.sbatch
#
# FULL RUN -- submit once PER COHORT (each is its own 1000-boot array):
#   sbatch --export=ALL,N_BOOT=10,COHORT=SPRINT run_autoc_bootstrap_array.sbatch
#   sbatch --export=ALL,N_BOOT=10,COHORT=ACCORD run_autoc_bootstrap_array.sbatch
#   sbatch --export=ALL,N_BOOT=10,COHORT=AASK   run_autoc_bootstrap_array.sbatch
#   sbatch --export=ALL,N_BOOT=10,COHORT=MDRD   run_autoc_bootstrap_array.sbatch
#   (COHORT=POOLED reproduces the pooled run.)
#
# Fragments are COHORT-tagged (frag_scalars_SPRINT_037.csv, ...), so all four
# can run at once without colliding. Combine each separately:
#   COHORT=SPRINT Rscript 03b_combine_bootstrap.R      (etc.)
#
# WATCH:   squeue -u $USER
# STATUS:  sacct -j <jobid> --format=JobID,State,ExitCode
# =============================================================================
