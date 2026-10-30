#!/bin/bash
#SBATCH --job-name=gtex_misid
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/logs/gtex_misid_%A_%a.out
#SBATCH --error=/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/logs/gtex_misid_%A_%a.err

# Submit with:  sbatch simulations/run_gtex.sh
# (the logs directory above must already exist, or the job fails silently)
# Afterwards:   Rscript simulations/Aggregate_gtex_results.R
#
# 240 scenarios split over 20 array tasks (12 each). Change --array and
# N_TASKS together; "%10" caps how many tasks run at once. A task that is
# requeued or times out resumes, because finished scenarios are skipped.

set -euo pipefail
module load r/4.5.1

cd /path/to/MISID_QC   # repo root on the cluster; the scripts use paths relative to it
export N_TASKS=20
# One thread per task: the parallelism is the array.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Extra runner arguments, e.g.
#   sbatch --export=ALL,RUN_ARGS="cv_folds=5 mismatch_repeats=1,3 target_specificity=0.9,0.95,0.99" simulations/run_gtex.sh
#   sbatch --export=ALL,RUN_ARGS="cv_folds=loo subset=rep==1" simulations/run_gtex.sh
# (no spaces or quotes inside RUN_ARGS values)
RUN_ARGS="${RUN_ARGS:-}"

# shellcheck disable=SC2086
Rscript simulations/Run_across_GTEX_simulations.R \
  task_id="${SLURM_ARRAY_TASK_ID}" n_tasks="${N_TASKS}" ${RUN_ARGS}
