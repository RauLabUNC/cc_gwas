#!/bin/bash

#SBATCH --job-name=test_with_tracking
#SBATCH --time=01:00:00
#SBATCH --mem=10000
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Set up logging
DATE_DIR=$(date +%Y-%m-%d)
RUN_ID=$(date +%H%M%S)_${SLURM_JOB_ID}_TRACK
LOG_DIR=".slurmlogs/${DATE_DIR}/${RUN_ID}"
mkdir -p "${LOG_DIR}"

# Redirect output
exec 1>"${LOG_DIR}/tracking.out" 2>"${LOG_DIR}/tracking.err"

# Load environment
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env
export LD_LIBRARY_PATH=/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib:$LD_LIBRARY_PATH

# Set test mode
export SNAKEMAKE_MODE="test"
export LOG_DIR="${LOG_DIR}"

# Parse arguments
PROFILE="${1:-test}"
DESCRIPTION="${2:-Test run with chr1 only and resource tracking}"

# Record start time
START_TIME=$(date +%Y-%m-%dT%H:%M:%S)
RUN_DATE=$(date +%Y%m%d_%H%M%S)

echo "================================================"
echo "Pipeline Run with Resource Tracking"
echo "================================================"
echo "Profile: $PROFILE"
echo "Start time: $START_TIME"
echo "Description: $DESCRIPTION"
echo "Log directory: $LOG_DIR"
echo "================================================"

# Use Snakemake directly
SNAKEMAKE_BIN="/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/bin/snakemake"

echo "Running Snakemake pipeline..."
${SNAKEMAKE_BIN} --snakefile scripts/from_scratch/snakemake/Snakefile \
    -j 10 \
    --cluster "sbatch --job-name smk_{rule} --time {resources.time} --mem={resources.mem_mb} -N 1 -n 1 -o /dev/null -e /dev/null --comment='{rule}:{wildcards}'" \
    --latency-wait 60 \
    --printshellcmds 2>&1 | tee "${LOG_DIR}/snakemake.log"

# Record end time
END_TIME=$(date +%Y-%m-%dT%H:%M:%S)

echo ""
echo "Pipeline completed at: $END_TIME"
echo "Collecting resource usage data..."

# Create report directory
REPORT_DIR="results/resource_reports/${RUN_DATE}_${PROFILE}"
mkdir -p "$REPORT_DIR"

# Wait for SLURM accounting database to update
echo "Waiting for SLURM accounting to update..."
sleep 30

# Collect SLURM usage data (using smk_ prefix now, not smk-)
./scripts/resource_tracking/collect_usage.sh \
    "$START_TIME" \
    "$END_TIME" \
    "smk_" \
    "${REPORT_DIR}/slurm_usage.tsv"

# Generate HTML report
echo "Generating resource usage report..."
Rscript scripts/resource_tracking/generate_report.R \
    "${REPORT_DIR}/slurm_usage.tsv" \
    "${REPORT_DIR}/usage_report.html"

# Copy benchmark files to report directory (they're small and useful)
echo "Copying benchmark data to report directory..."
if [ -d "bench" ]; then
    cp -r bench "${REPORT_DIR}/" 2>/dev/null || true
fi

# Create a symlink to the full logs instead of copying
ln -s "$(pwd)/${LOG_DIR}" "${REPORT_DIR}/full_logs" 2>/dev/null || true

echo ""
echo "================================================"
echo "Resource tracking complete!"
echo "Report generated at: ${REPORT_DIR}/usage_report.html"
echo "================================================"