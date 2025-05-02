#!/bin/bash

#SBATCH --job-name=snake_norm
#SBATCH --time=24:00:00   # Reduced time, adjust as needed for R script
#SBATCH --mem=4000      # Adjusted memory, maybe less needed? Check usage.
#SBATCH --output=.slurmlogs/norm_%j.out # Separate log directory or prefix

# Load necessary modules
# Ensure R and necessary packages (tidyverse, MASS) are available
module load r r/4.3.1 # Or your specific R version/module


# Run Snakemake using the specific Snakefile_boxCox
# Executes the 'all' rule by default
snakemake --snakefile Snakefile_boxCox \
          -j 250 \
          --latency-wait 60 \
          --cluster "sbatch --time 16:00:00 --mem=4000 -N 1 -n 1 -o .slurmlogs/norm_job_%j.out" `# Cluster settings for individual R script jobs` \
          --keep-going \
          --rerun-incomplete 