#!/bin/bash
#SBATCH --job-name=make_counts
#SBATCH --time=00:30:00
#SBATCH --mem=16000
#SBATCH --output=.slurmlogs/misc/%j.out
#SBATCH --error=.slurmlogs/misc/%j.err

set -euo pipefail

# Load conda exactly like the Snakemake submission scripts
set +u
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env
set -u

# Match environment library handling from working pipeline scripts
export LD_LIBRARY_PATH=/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib:$LD_LIBRARY_PATH
unset R_HOME || true

echo "===================================="
echo "miQTL Environment Loaded"
echo "===================================="
echo "Using Rscript: $(which Rscript)"
Rscript -e 'cat("R version:", R.version.string, "\n"); cat("LibPaths:\n"); print(.libPaths())'

# Run script with the env Rscript
Rscript scripts/03_makeCounts.R