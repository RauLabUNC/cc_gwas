#!/bin/bash
# Script to rebuild the miqtl environment from scratch

set -euo pipefail

echo "Rebuilding miqtl-env environment..."

# Load mamba
source ~/mambaforge/etc/profile.d/conda.sh

# Remove old environment if exists
mamba env remove -n miqtl-env -y || true

# Create new environment from YAML
mamba env create -f envs/environment.yml

# Activate environment
conda activate miqtl-env

# Install miqtl from GitHub
echo "Installing miqtl from GitHub..."
Rscript -e "if (!require('remotes')) install.packages('remotes'); remotes::install_github('gkeele/miqtl')"

# Export frozen environment
echo "Saving environment snapshot..."
conda env export --no-builds > envs/environment-frozen.yml
mamba list --explicit > envs/conda-spec.txt

echo "Environment rebuilt successfully!"
echo "Activate with: conda activate miqtl-env"