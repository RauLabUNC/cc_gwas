#!/bin/bash
# create_comparison_report.sh - Create a resource comparison report from job logs

set -euo pipefail

LOG_DIR="${1:-.slurmlogs/2025-09-02/145732_14065877_TRACK}"
OUTPUT="${2:-resource_comparison.txt}"

echo "==================================================================" > "$OUTPUT"
echo "Resource Usage Comparison Report" >> "$OUTPUT"
echo "==================================================================" >> "$OUTPUT"
echo "" >> "$OUTPUT"

# Find all time files
for time_file in $(find "$LOG_DIR" -name "time" -type f); do
    # Extract job info from path
    job_path=$(dirname "$time_file")
    rule=$(basename $(dirname "$job_path"))
    wildcards=$(basename "$job_path")
    
    echo "Rule: $rule" >> "$OUTPUT"
    echo "Job: $wildcards" >> "$OUTPUT"
    echo "------------------------------------------------------------------" >> "$OUTPUT"
    
    # Get SLURM job ID from time file
    job_id=$(grep "SLURM_JOB_ID=" "$time_file" 2>/dev/null | cut -d= -f2)
    
    # Get actual usage from time file
    if [ -f "$time_file" ]; then
        cpu_time=$(grep "User time" "$time_file" | awk '{print $4}')
        wall_time=$(grep "Elapsed (wall clock)" "$time_file" | awk '{print $8}')
        max_rss=$(grep "Maximum resident" "$time_file" | awk '{print $6}')
        cpu_percent=$(grep "Percent of CPU" "$time_file" | awk '{print $7}')
        
        echo "ACTUAL USAGE:" >> "$OUTPUT"
        echo "  Wall time: $wall_time" >> "$OUTPUT"
        echo "  CPU time: ${cpu_time}s" >> "$OUTPUT"
        echo "  Max memory: $((max_rss/1024)) MB" >> "$OUTPUT"
        echo "  CPU utilization: $cpu_percent" >> "$OUTPUT"
    fi
    
    # Get requested resources from SLURM if job ID exists
    if [ ! -z "$job_id" ] && [ "$job_id" != "" ]; then
        echo "" >> "$OUTPUT"
        echo "REQUESTED (Job $job_id):" >> "$OUTPUT"
        sacct -j "$job_id" --format=Timelimit,ReqMem,AllocCPUS -n | head -1 | \
            awk '{print "  Time limit: " $1 "\n  Memory: " $2 "\n  CPUs: " $3}' >> "$OUTPUT"
        
        # Calculate efficiency
        echo "" >> "$OUTPUT"
        echo "EFFICIENCY:" >> "$OUTPUT"
        
        # Get SLURM's view
        sacct -j "$job_id.batch" --format=MaxRSS,CPUTimeRAW,Elapsed -n | head -1 | \
            awk '{
                if ($1 != "") {
                    rss_mb = $1
                    gsub(/K/, "", rss_mb)
                    rss_mb = rss_mb / 1024
                    print "  Memory used: " rss_mb " MB (SLURM measurement)"
                }
                print "  CPU time (SLURM): " $2 " seconds"
                print "  Actual runtime: " $3
            }' >> "$OUTPUT"
    fi
    
    echo "" >> "$OUTPUT"
    echo "==================================================================" >> "$OUTPUT"
    echo "" >> "$OUTPUT"
done

# Also check benchmark files if they exist
if [ -d "bench" ]; then
    echo "SNAKEMAKE BENCHMARK DATA:" >> "$OUTPUT"
    echo "------------------------------------------------------------------" >> "$OUTPUT"
    for bench_file in $(find bench -name "*.txt" -type f); do
        echo "File: $bench_file" >> "$OUTPUT"
        cat "$bench_file" >> "$OUTPUT"
        echo "" >> "$OUTPUT"
    done
fi

echo "Report saved to: $OUTPUT"
cat "$OUTPUT"