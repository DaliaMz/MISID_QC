#!/bin/bash
#SBATCH --job-name=sims_misid
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/logs/sims_misid_%j.out
#SBATCH --error=/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/logs/sims_misid_%j.err

# Submit with:  sbatch simulations/run_sims.sh
# (the logs directory above must already exist, or the job fails silently)
#
# One job, parallel inside R (mclapply forks, so the simulated datasets are
# loaded once and shared). It is a single job rather than an array because
# the runner writes one summary file at the end. Finished scenarios are
# skipped, so if the job times out just resubmit it.
# Keep --cpus-per-task and cores in step.

set -euo pipefail
module load r/4.5.1

cd /scratch/project_mnt/S0138/DMizikovsky/misid/MISID_QC   # repo root on the cluster; the scripts use paths relative to it
CORES="${SLURM_CPUS_PER_TASK:-8}"
# One BLAS thread per worker: the parallelism is the forked workers.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Extra runner arguments, e.g.
#   sbatch --export=ALL,RUN_ARGS="cv_folds=5 mismatch_repeats=1,3 target_specificity=0.9,0.95,0.99" simulations/run_sims.sh
#   sbatch --export=ALL,RUN_ARGS="cv_folds=loo subset=rep==1" simulations/run_sims.sh
# (no spaces or quotes inside RUN_ARGS values)
RUN_ARGS="${RUN_ARGS:-}"

# shellcheck disable=SC2086
Rscript simulations/Run_across_simulated_data.R cores="${CORES}" ${RUN_ARGS}
