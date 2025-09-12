# Genetic and transcriptomic determinants of heart failure susceptibility and response  
This doc describes a bioinformatics workflow that joins genetic, transcriptomic, and phenotypic measures from mouse cardiomyopathy models. As (hopefully!) described in Kimball *et al* 2026, we induced heart failure in 63 diverse strains of mice, reprenting 400+ mice total. We did this to learn, for instance, how there are heritable variations of the DNA which seem to minimize the how much the heart remodels after insult, minimizing stress later. 

### What we did:


### Full Pipeline
```bash
# Run full pipeline
sbatch scripts//snakemake/.sh
```
## Pipeline Structure
The scripts are numbered to show execution order:
```
scripts/from_scratch/
├── snakemake/
│   ├── 01_ropScan.R         # Genome-wide QTL scan
│   ├── 02_permTest.R         # Permutation thresholds
│   └── Snakefile             # Main workflow (complete pipeline)
├── 03_detectSigLoci.R        # Find significant regions
├── 04_joinLoci.R             # Merge trait & expression QTLs (miQTL only)
├── 10_joinPositionID.R       # Consolidate unique positions
├── 11_mergeLociDetectGenes.R # Find genes in loci
├── 12_splitGenesToSpecies.R  # Split mouse/human genes
├── 13_multTrait_cis-eQTL_nrvmExp.R # Gene list (modified - no PyLMM)
├── 20_getPathogenicity.R     # Disease associations from Open Targets
├── 21_getMouseGenePhenotypes.R # Mouse phenotypes (BLOCKED - MouseMine offline)
├── 22_makeLociPackets.R     # Final QTL packets with plots (no PyLMM)
└── notes.md                  # Detailed docs for each script
```
Logs are organized in `.slurmlogs/YYYY-MM-DD/HHMMSS_JOBID/` with separate files for each rule.

## Environment Management

```bash
# Activate the consolidated environment
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env

# Rebuild from scratch if needed
./scripts/rebuild_env.sh

# Environment specification in envs/environment.yml
```
The miqtl-env contains R 4.4.0, Python 3.11, Snakemake 7.32.4, and all necessary R packages. 