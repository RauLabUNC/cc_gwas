#!/bin/bash

#SBATCH --job-name=detect_loci_only
#SBATCH --time=01:00:00      
#SBATCH --mem=20000           

# Set up organized logging structure FIRST, before any output
DATE_DIR=$(date +%Y-%m-%d)
RUN_ID=$(date +%H%M%S)_${SLURM_JOB_ID}
LOG_DIR=".slurmlogs/${DATE_DIR}/${RUN_ID}"
mkdir -p "${LOG_DIR}/jobs"

# Redirect all output to date-specific directory
exec 1>"${LOG_DIR}/main.out" 2>"${LOG_DIR}/main.err"

# Load conda/mamba for miqtl environment (contains everything we need)
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env

# Add library path for ICU libraries needed by RCurl/InterMineR
export LD_LIBRARY_PATH=/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib:$LD_LIBRARY_PATH

echo "================================================"
echo "QTL Pipeline Run: ${RUN_ID}"
echo "================================================"
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Log directory: ${LOG_DIR}"
echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo "R: $(which R)"
echo "Snakemake: $(which snakemake)"

# Use the conda environment's Snakemake explicitly
SNAKEMAKE_BIN="/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/bin/snakemake"
echo "Using Snakemake: ${SNAKEMAKE_BIN}"
echo "Snakemake version: $(${SNAKEMAKE_BIN} --version)"
echo "================================================"
echo ""

# Run the pipeline with organized logging and named jobs
${SNAKEMAKE_BIN} --snakefile scripts/from_scratch/snakemake/Snakefile \
          -j 10 \
          --cluster "sbatch --job-name {rule} --time {resources.time} --mem={resources.mem_mb} -N 1 -n 1 -o ${LOG_DIR}/jobs/{rule}_%j.out -e ${LOG_DIR}/jobs/{rule}_%j.err" \
          --latency-wait 60 \
          --printshellcmds 2>&1 | tee "${LOG_DIR}/snakemake.log"

echo ""
echo "Pipeline completed at: $(date)"