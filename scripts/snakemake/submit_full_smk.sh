#!/bin/bash
#SBATCH --job-name=full_permutation_run
#SBATCH --time=012:00:00
#SBATCH --mem=4000
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

DATE_DIR=$(date +%Y-%m-%d)
RUN_ID="$(date +%H%M%S)_${SLURM_JOB_ID:-parent}"
LOG_DIR=".slurmlogs/${DATE_DIR}/${RUN_ID}"
mkdir -p "${LOG_DIR}/jobs" "${LOG_DIR}/bench"

SACCT_START_ISO=$(date +%Y-%m-%dT%H:%M:%S)

exec 1>"${LOG_DIR}/main.out" 2>"${LOG_DIR}/main.err"

# Conda hooks don’t like `set -u`
set +u
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env
set -u

export LD_LIBRARY_PATH=/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib:$LD_LIBRARY_PATH
export SNAKEMAKE_MODE="full"
export RUN_ID LOG_DIR

echo "================================================"
echo "QTL Pipeline FULL Run: ${RUN_ID}"
echo "================================================"
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Log directory: ${LOG_DIR}"
echo "Python: $(which python)"; python --version
echo "R: $(which R)"
SNAKEMAKE_BIN="/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/bin/snakemake"
echo "Using Snakemake: ${SNAKEMAKE_BIN}"
echo "Snakemake version: $(${SNAKEMAKE_BIN} --version)"
echo "================================================"
echo ""

${SNAKEMAKE_BIN} --snakefile scripts/snakemake/smk_para_perm.smk \
  -j 500 \
  --rerun-incomplete --keep-going \
  --latency-wait 60 \
  --cluster "sbatch \
    --job-name=${RUN_ID}_{rule}-{wildcards} --comment={rule}:{wildcards} \
    --time={resources.time} \
    --mem={resources.mem_mb} \
    -N 1 -n 1 --cpus-per-task={threads} \
    -o ${LOG_DIR}/jobs/{rule}/shell_logs/%j.out \
    -e ${LOG_DIR}/jobs/{rule}/shell_logs/%j.err" \
  --envvars RUN_ID LOG_DIR SNAKEMAKE_MODE \
  --stats  "${LOG_DIR}/snakemake_stats.json" \
  --printshellcmds 2>&1 | tee "${LOG_DIR}/snakemake.log"
  # --report "${LOG_DIR}/snakemake_report.html"                              # left commented

SNKM_STATUS=${PIPESTATUS[0]}

# NEW: capture an explicit end time and give accounting a brief moment
SACCT_END_ISO=$(date +%Y-%m-%dT%H:%M:%S)
sleep 120

# NEW: let sacct filter by name; widen JobName to avoid truncation surprises
sacct -X -P -n -S "${SACCT_START_ISO}" -E "${SACCT_END_ISO}" \
  --name="${RUN_ID}_%" \
  -o JobID,JobName%60,Partition,AllocCPUS,Elapsed,State,AveCPU,CPUTimeRAW,ReqMem,MaxRSS,MaxVMSize,ExitCode \
  > "${LOG_DIR}/slurm_metrics.tsv" || true

echo ""
echo "Pipeline completed at: $(date)  (status=${SNKM_STATUS})"
exit "${SNKM_STATUS}"

