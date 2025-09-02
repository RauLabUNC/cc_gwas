#!/usr/bin/env bash
# run_with_tracking.sh - Run Snakemake pipeline with automatic resource tracking

set -euo pipefail

# Parse arguments
PROFILE="${1:-test}"  # "test" or "full"
DESCRIPTION="${2:-Pipeline run}"

# Set submission script based on profile
if [ "$PROFILE" = "test" ]; then
    SUBMIT_SCRIPT="scripts/from_scratch/snakemake/submit_test_mode.sh"
    echo "Running TEST pipeline with resource tracking"
else
    SUBMIT_SCRIPT="scripts/from_scratch/snakemake/submit_with_existing.sh"
    echo "Running FULL pipeline with resource tracking"
fi

# Record start time
START_TIME=$(date +%Y-%m-%dT%H:%M:%S)
RUN_DATE=$(date +%Y%m%d_%H%M%S)

echo "================================================"
echo "Pipeline Run with Resource Tracking"
echo "================================================"
echo "Profile: $PROFILE"
echo "Start time: $START_TIME"
echo "Description: $DESCRIPTION"
echo "================================================"

# Submit the pipeline job and capture the job ID
echo "Submitting pipeline..."
JOB_OUTPUT=$(sbatch --parsable $SUBMIT_SCRIPT)
MAIN_JOB_ID=$(echo $JOB_OUTPUT | awk '{print $1}')

echo "Main job ID: $MAIN_JOB_ID"
echo "Waiting for pipeline to complete..."

# Wait for the main job to finish
while squeue -j $MAIN_JOB_ID 2>/dev/null | grep -q $MAIN_JOB_ID; do
    sleep 30
done

# Record end time
END_TIME=$(date +%Y-%m-%dT%H:%M:%S)

echo "Pipeline completed at: $END_TIME"
echo ""
echo "Collecting resource usage data..."

# Create report directory
REPORT_DIR="results/resource_reports/${RUN_DATE}_${PROFILE}"
mkdir -p "$REPORT_DIR"

# Collect SLURM usage data
./scripts/resource_tracking/collect_usage.sh \
    "$START_TIME" \
    "$END_TIME" \
    "smk-" \
    "${REPORT_DIR}/slurm_usage.tsv"

# Generate HTML report
echo "Generating resource usage report..."
Rscript scripts/resource_tracking/generate_report.R \
    "${REPORT_DIR}/slurm_usage.tsv" \
    "${REPORT_DIR}/usage_report.html"

# Copy logs and benchmarks to report directory
echo "Collecting logs and benchmarks..."
if [ -d "logs" ]; then
    cp -r logs "${REPORT_DIR}/" 2>/dev/null || true
fi
if [ -d "bench" ]; then
    cp -r bench "${REPORT_DIR}/" 2>/dev/null || true
fi

# Create summary file
cat > "${REPORT_DIR}/run_summary.txt" << EOF
Pipeline Run Summary
====================
Profile: $PROFILE
Description: $DESCRIPTION
Start Time: $START_TIME
End Time: $END_TIME
Main Job ID: $MAIN_JOB_ID
Report Directory: $REPORT_DIR

Files Generated:
- usage_report.html: Interactive HTML report with resource usage analysis
- slurm_usage.tsv: Raw SLURM accounting data
- logs/: Individual job logs with time measurements
- bench/: Snakemake benchmark files
- plots/: Resource usage visualizations

To view the report:
  firefox ${REPORT_DIR}/usage_report.html
EOF

echo ""
echo "================================================"
echo "Resource tracking complete!"
echo "Report generated at: ${REPORT_DIR}/usage_report.html"
echo "================================================"
echo ""
echo "Quick summary:"
sacct -j $MAIN_JOB_ID --format=JobID,JobName%30,State,Elapsed,MaxRSS,ExitCode

# Open report if display is available
if [ -n "${DISPLAY:-}" ]; then
    echo ""
    echo "Opening report in browser..."
    xdg-open "${REPORT_DIR}/usage_report.html" 2>/dev/null || \
    firefox "${REPORT_DIR}/usage_report.html" 2>/dev/null || \
    echo "Please open ${REPORT_DIR}/usage_report.html in your browser"
fi