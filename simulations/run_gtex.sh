#!/bin/bash
#SBATCH --job-name=gtex_misid
#SBATCH --array=1-30%15
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/home/projects/ssi_gen1/people/dalmiz/Sample_Misassignment_0826/MISID_QC/simulations/logs/gtex_misid_%A_%a.out
#SBATCH --error=/home/projects/ssi_gen1/people/dalmiz/Sample_Misassignment_0826/MISID_QC/simulations/logs/gtex_misid_%A_%a.err

# Submit with:  sbatch simulations/run_gtex.sh
# (simulations/logs must already exist, or the job fails silently)
# Afterwards:   Rscript simulations/Aggregate_gtex_results.R
#
# Up to 300 scenarios (150 base settings with and without unknown classes,
# random tissue sets, mislabel rates 0 to 0.4) split over 30 array tasks
# (10 each). Change --array and
# N_TASKS together; "%15" caps how many tasks run at once. A task that is
# requeued or times out resumes, because finished scenarios are skipped.
# The default run is the primary config only (5-fold CV, 3 label-swap repeats,
# margin_coverage 0.025, target_specificity 0.95). Run one-at-a-time sweeps on a
# subset afterwards, e.g. RUN_ARGS="cv_folds=5 mismatch_repeats=1,3,5 subset=rep==1".
# Setting RUN_ARGS replaces the default entirely. The agreement gate is off
# (min_agreement=0); first_pass_agreement is recorded for each scenario. Finished
# configs are skipped, so use a new out= directory or force=TRUE after a code change.

set -euo pipefail
module load r/4.5.1

cd /home/projects/ssi_gen1/people/dalmiz/Sample_Misassignment_0826/MISID_QC   # repo root on the cluster; the scripts use paths relative to it
export N_TASKS=30
# One thread per task: the parallelism is the array.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Extra runner arguments, e.g.
#   sbatch --export=ALL,RUN_ARGS="cv_folds=5 mismatch_repeats=3 margin_coverage=0.01,0.025,0.05 subset=rep==1" simulations/run_gtex.sh
#   sbatch --export=ALL,RUN_ARGS="cv_folds=loo subset=rep==1" simulations/run_gtex.sh
# (no spaces or quotes inside RUN_ARGS values)
RUN_ARGS="${RUN_ARGS:-cv_folds=5}"

# shellcheck disable=SC2086
Rscript simulations/Run_across_GTEX_simulations.R \
  task_id="${SLURM_ARRAY_TASK_ID}" n_tasks="${N_TASKS}" ${RUN_ARGS}
