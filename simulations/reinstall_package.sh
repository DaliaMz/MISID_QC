#!/bin/bash
set -euo pipefail

export HOME=/home/projects/ssi_gen1/people/dalmiz
source ~/miniforge3/bin/activate
conda activate jupyter

cd "$(dirname "${BASH_SOURCE[0]}")/.."
R CMD INSTALL .
Rscript -e 'cat("Loaded package from:", find.package("SampleAnnotationLikelihood"), "\n")'
