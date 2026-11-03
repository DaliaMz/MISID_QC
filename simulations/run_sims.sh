#!/bin/bash
#SBATCH --job-name=sims_misid
#SBATCH --array=1-60%25
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/logs/sims_misid_%A_%a.out
#SBATCH --error=/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/logs/sims_misid_%A_%a.err

# Submit with:  sbatch simulations/run_sims.sh
# (the logs directory above must already exist, or the job fails silently)
# Afterwards:   Rscript simulations/Run_across_simulated_data.R combine=TRUE
#               (writes simulations/results/summary_v2.rds and .csv)
#
# The 360 simulated scenarios are split over 60 array tasks (task k takes
# scenarios k, k+60, ..., interleaved so slow scenarios spread evenly). Change
# --array and N_TASKS together; "%25" caps how many tasks run at once. Each task
# writes its own row files, so tasks never write to the same file, and a task
# that is requeued or times out resumes because finished scenarios are skipped.
# Each task reads the whole simulated dataset file once, hence the memory request.
#
# The default run is the primary config only (5-fold CV, 3 label-swap repeats,
# target_specificity 0.95) with the agreement gate off (min_agreement=0);
# first_pass_agreement is recorded per scenario. Finished configs are skipped,
# so use a new out= directory or force=TRUE after a code change. Setting RUN_ARGS
# replaces the default arguments entirely.

set -euo pipefail
module load r/4.5.1

cd /scratch/project_mnt/S0138/DMizikovsky/misid/MISID_QC   # repo root on the cluster; the scripts use paths relative to it
export N_TASKS=60
# One thread per task: the parallelism is the array.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Extra runner arguments, e.g. one-at-a-time sweeps on a subset
#   sbatch --export=ALL,RUN_ARGS="cv_folds=5 mismatch_repeats=1,3,5 subset=rep==1" simulations/run_sims.sh
#   sbatch --export=ALL,RUN_ARGS="cv_folds=5 target_specificity=0.9,0.95,0.99 subset=rep==1" simulations/run_sims.sh
#   sbatch --export=ALL,RUN_ARGS="cv_folds=loo subset=rep==1" simulations/run_sims.sh
# (no spaces or quotes inside RUN_ARGS values)
RUN_ARGS="${RUN_ARGS:-cv_folds=5}"

# shellcheck disable=SC2086
Rscript simulations/Run_across_simulated_data.R \
  task_id="${SLURM_ARRAY_TASK_ID}" n_tasks="${N_TASKS}" ${RUN_ARGS}
