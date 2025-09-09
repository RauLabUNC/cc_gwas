#!/usr/bin/env bash
# generate_report_from_logs.sh - Generate resource report from existing SLURM logs

set -euo pipefail

# Usage message
if [ $# -lt 2 ]; then
    echo "Usage: $0 <start_time> <end_time> [job_prefix] [output_dir]"
    echo "Example: $0 '2025-09-02T14:00:00' '2025-09-02T15:00:00'"
    echo "Example: $0 '2 hours ago' 'now' smk- my_report"
    exit 1
fi

# Parse arguments
START_TIME="$1"
END_TIME="$2"
JOB_PREFIX="${3:-smk-}"
OUTPUT_DIR="${4:-results/resource_reports/$(date +%Y%m%d_%H%M%S)}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "================================================"
echo "Generating Resource Report from SLURM Logs"
echo "================================================"
echo "Time range: $START_TIME to $END_TIME"
echo "Job prefix: $JOB_PREFIX"
echo "Output directory: $OUTPUT_DIR"
echo "================================================"

# Collect SLURM data
echo "Collecting SLURM accounting data..."
./scripts/resource_tracking/collect_usage.sh \
    "$START_TIME" \
    "$END_TIME" \
    "$JOB_PREFIX" \
    "${OUTPUT_DIR}/slurm_usage.tsv"

# Check if we got any data
JOB_COUNT=$(tail -n +2 "${OUTPUT_DIR}/slurm_usage.tsv" | wc -l)
if [ "$JOB_COUNT" -eq 0 ]; then
    echo "WARNING: No jobs found matching criteria"
    echo "Check your time range and job prefix"
    exit 1
fi

echo "Found $JOB_COUNT jobs"

# Generate HTML report
echo "Generating HTML report..."
Rscript scripts/resource_tracking/generate_report.R \
    "${OUTPUT_DIR}/slurm_usage.tsv" \
    "${OUTPUT_DIR}/usage_report.html"

echo ""
echo "================================================"
echo "Report generated successfully!"
echo "Location: ${OUTPUT_DIR}/usage_report.html"
echo "================================================"

# Show quick summary
echo ""
echo "Quick resource summary:"
echo "------------------------"
tail -n +2 "${OUTPUT_DIR}/slurm_usage.tsv" | \
    awk -F'|' '{
        split($6, elapsed, ":");
        total_time += elapsed[1]*3600 + elapsed[2]*60 + elapsed[3];
        n++;
    } END {
        if (n > 0) {
            hours = int(total_time / 3600);
            minutes = int((total_time % 3600) / 60);
            printf "Total jobs: %d\n", n;
            printf "Total wallclock time: %02d:%02d hours\n", hours, minutes;
        }
    }'