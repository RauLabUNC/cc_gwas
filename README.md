# CC miQTL Pipeline Rebuild

This is the rebuilt version of the Collaborative Cross miQTL analysis pipeline. The goal is to get everything from the original scattered scripts into a single, reproducible workflow.

## Current Status (September 2, 2025)

### Latest Updates
- ✅ **Test mode implemented**: 5-minute pipeline validation (chr1 only, 5 permutations)
- ✅ **Resource tracking system**: Complete monitoring with HTML reports and SLURM integration
- ✅ **Optimized resource allocation**: Reduced memory from 6GB to 1GB for test jobs
- ✅ **Git repository initialized**: Connected to GitHub at https://github.com/RauLabUNC/cc_gwas
- 🔄 **In progress**: Removing hardcoded paths and integrating data generation

## Quick Start

### Test Mode (5 minutes)
```bash
# Run test with resource tracking
sbatch scripts/resource_tracking/submit_with_tracking.sh test "My test description"

# Or without tracking
sbatch scripts/from_scratch/snakemake/submit_test_mode.sh
```

### Full Pipeline
```bash
# Run full pipeline
sbatch scripts/from_scratch/snakemake/submit_with_existing.sh
```

### Resource Reports
After runs complete, find reports in:
- `results/resource_reports/[timestamp]_test/usage_report.html`
- `results/resource_reports/[timestamp]_test/comparison.txt`

## Previous Status (Aug 29, 2025)

**MAJOR PROGRESS**: Complete pipeline configured from QTL scanning through final loci packet generation!

Key achievements today:
- Removed all PyLMM dependencies (focusing on miQTL results only)
- Integrated makeLociPackets.R as final pipeline step
- Fixed circular dependencies in Snakefile
- Symlinked bulk expression and NRVM data from original project
- Fixed SLURM logging to only save in date-specific directories
- Added LD_LIBRARY_PATH workaround for InterMineR

**BLOCKING ISSUE**: MouseMine appears to be offline/retired - InterMineR cannot connect to it for phenotype annotations (script 21).

## What's Working

- **Genome scans**: The ROP scan from miQTL package works with the symlinked genome cache
- **Permutation testing**: Cut down to 20 permutations for testing (normally 50-100)
- **Loci detection**: Identifies significant regions based on LOD thresholds
- **Gene annotation**: Pulls genes within QTL regions and gets disease associations from Open Targets

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

## Quick Start

```bash
# Activate environment
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env

# Run with existing scan results (recommended - saves hours)
sbatch scripts/from_scratch/snakemake/submit_with_existing.sh

# Or run full pipeline from scratch
sbatch scripts/from_scratch/snakemake/submit_test_chr1_v2.sh
```

Logs are organized in `.slurmlogs/YYYY-MM-DD/HHMMSS_JOBID/` with separate files for each rule.

## Data Paths

- **Phenotypes**: `data/processed/phenotypes/boxcox_individual_{Ctrl,Iso}.csv`
- **Genome cache**: Symlinked to `/proj/raulab/projects/cc_gwas/data/raw/genomes/`
- **Test outputs**: `data/processed/test_chr1/`

## Next Steps

1. ✅ ~~Install snakemake slurm executor~~ → Using Snakemake 7.32.4 instead
2. ✅ ~~Add remaining scripts to workflow~~ → All pipeline steps integrated
3. Fix InterMineR/RCurl library dependencies for script 21
4. Eventually integrate the NRVM expression data and cis-eQTL analysis

## Known Issues

- **MouseMine offline**: InterMineR cannot connect to MouseMine database (appears to be retired/unavailable as of Aug 29, 2025)
- **LD_LIBRARY_PATH hardcoded**: ICU libraries path hardcoded in submission script - needs dynamic solution
- **PyLMM removed**: All PyLMM dependencies removed from pipeline - using miQTL results only
- Scripts process all available scan files, not just the test traits

## Original Pipeline

The original messy version is in `/proj/raulab/projects/cc_gwas/` if you need to reference something. Main difference is we're trying to get everything into a single workflow instead of manually running 15 different scripts in order.

## Environment Management

```bash
# Activate the consolidated environment
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env

# Rebuild from scratch if needed
./scripts/rebuild_env.sh

# Environment specification in envs/environment.yml
```

The miqtl-env contains R 4.4.0, Python 3.11, Snakemake 7.32.4, and all necessary R packages. The miqtl package is installed from GitHub to `/proj/raulab/users/brian/cc_gwas/R_packages/`.