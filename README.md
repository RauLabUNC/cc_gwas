# Genetic and transcriptomic determinants of heart failure susceptibility and response  
This doc describes a bioinformatics workflow that joins genetic, transcriptomic, and phenotypic measures from mouse cardiomyopathy models. As (hopefully!) described in Kimball *et al* 2026, we induced heart failure in 63 diverse strains of mice, reprenting 400+ mice total. We did this to learn, for instance, how there are heritable variations of the DNA which seem to minimize the how much the heart remodels after insult, minimizing stress later. 

### What we did:


### Full Pipeline
```bash
# Run full pipeline
sbatch scripts/snakemake/.sh
```
## Pipeline Structure

Split into four parts:
- `scripts/align` STAR alignment of bulk RNAseq
- `scripts/traitQTL` haplotype-trait associations with QTL mapping
- `scripts/eQTL` haploytype-expression QTL mapping 
- `scripts/packets ` Join all results into summary packets

Descriptions on how to produce each are found within their respective directories.


Logs are organized in `.slurmlogs/YYYY-MM-DD/HHMMSS_JOBID/` with separate files for each rule.

## Environment Management

```bash
# for troubleshooting
srun -t 5:00:00 -p interact -n 1 --cpus-per-task=1 --mem=16g  --pty /bin/bash
# Activate the consolidated environment
set +u
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env
set -u

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
R


# Rebuild from scratch if needed
./scripts/rebuild_env.sh

# Environment specification in envs/environment.yml
```
The miqtl-env contains R 4.4.0, Python 3.11, Snakemake 7.32.4, and all necessary R packages. 