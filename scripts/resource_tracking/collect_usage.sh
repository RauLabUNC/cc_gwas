#!/usr/bin/env bash
# collect_usage.sh - Collect SLURM job resource usage data
# Usage:
#   ./collect_usage.sh [START_ISO] [END_ISO] [JOBNAME_PREFIX] [OUT_TSV]
# Examples:
#   ./collect_usage.sh                         # today 00:00 → now, prefix smk-
#   ./collect_usage.sh 2025-09-02T08:00:00     # explicit start time

set -euo pipefail

# Parse arguments with defaults
START_ISO="${1:-$(date -d 'today 00:00:00' +%Y-%m-%dT%H:%M:%S 2>/dev/null || date -v0H -v0M -v0S +%Y-%m-%dT%H:%M:%S)}"
END_ISO="${2:-$(date +%Y-%m-%dT%H:%M:%S)}"
PREFIX="${3:-smk-}"
OUT="${4:-results/resource_reports/slurm_usage.tsv}"

# Create output directory if needed
mkdir -p "$(dirname "$OUT")"

# Define header for the output TSV
HEADER="JobID|JobName|State|Start|End|Elapsed|AllocCPUS|ReqCPUS|ReqMem|Timelimit|CPUTimeRAW|UserCPU|SystemCPU|AveCPU|MaxRSS|AveRSS|AveVMSize|AveDiskRead|AveDiskWrite|ExitCode"
echo "$HEADER" > "$OUT"

# Collect data from SLURM accounting
sacct -X -P -n -S "$START_ISO" -E "$END_ISO" --name="${PREFIX}%" \
  --format=JobID,JobName,State,Start,End,Elapsed,AllocCPUS,ReqCPUS,ReqMem,Timelimit,CPUTimeRAW,UserCPU,SystemCPU,AveCPU,MaxRSS,AveRSS,AveVMSize,AveDiskRead,AveDiskWrite,ExitCode \
  >> "$OUT"

echo "Wrote $OUT (window: $START_ISO .. $END_ISO; prefix: $PREFIX)" >&2
echo "Found $(wc -l < "$OUT") jobs (including header)" >&2