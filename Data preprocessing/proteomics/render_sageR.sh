#!/bin/bash
#SBATCH -J rmd_render
#SBATCH -p node2
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 04:00:00
#SBATCH --mem=8G
#SBATCH -o rmd.%j.out
#SBATCH -e rmd.%j.err

set -euo pipefail

# --- Paths ---
WORKDIR="/data/Storage_DAS03/home3/menteeF/yeongyu/Proteomics"
RMD="sageR_processing.Rmd"

# --- Environment ---
module load R/4.2.1
export R_LIBS_USER="$HOME/R/4.2/library"
mkdir -p "$R_LIBS_USER"

# Avoid runaway threading from BLAS/OpenMP during knit
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "JobID=${SLURM_JOB_ID:-NA} Host=$(hostname) PWD=$WORKDIR"
echo "Rscript: $(which Rscript)"; Rscript --version

cd "$WORKDIR"

# Knit to HTML (PDF needs LaTeX, so keep it HTML unless you’ve set that up)
Rscript -e 'rmarkdown::render("'"$RMD"'", output_format="html_document")'
echo "[OK] render finished: ${RMD/.Rmd/.html}"
