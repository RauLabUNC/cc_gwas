# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics project for miQTL (micro-RNA Quantitative Trait Loci) analysis using the Collaborative Cross (CC) mouse model. The project analyzes gene expression data to identify genetic loci associated with expression traits under control and isoproterenol (iso) treatment conditions.

## Current State & Evolution

### Key Scripts Overview

1. **Everything.R** - A consolidated pipeline script that combines multiple analysis steps:
   - Expression data processing (Steps 0-5): Sample filtering, batch correction (ComBat), CPM filtering, VST transformation
   - ROP QTL scanning 
   - Permutation threshold determination
   - Significant haplotype block identification
   - Data consolidation across all genes
   - cis-eQTL analysis within trait loci

2. **Brian_loci_packets.R** - The original comprehensive loci analysis pipeline that generates "QTL packets":
   - Loads extensive relational data (mouse/human phenotypes, orthology, disease associations)
   - Groups overlapping loci into clusters using graph theory
   - Generates locus zoom plots with multiple tracks (miQTL, pyLMM, founder effects, genes)
   - Creates detailed Excel reports with gene information, cardiac trait associations
   - Produces README summaries focused on cardiac-related genes
   - Output: `/proj/raulab/users/Anh/CC-miQTL-clean/Results/QTLpackets/`

3. **Anh_loci_packets.R** - Modified version by Anh with key changes:
   - Different path for merged gene info: `/proj/raulab/users/Anh/CC-miQTL-clean/Results/QTLpackets/multTrait_cis-eQTL_nrvmExp_250523Anh.csv`
   - Split heart pattern into separate human and mouse patterns for more precise matching
   - Otherwise maintains same structure as Brian's version

4. **Update_loci_tabs.R** - Post-processing script to enhance existing QTL packets:
   - Adds CPM values, log2 fold changes, DESeq2 p-values
   - Adds cis-eQTL flags for control and iso conditions
   - Updates existing Excel files in place

### Key Changes Since Previous Work

- The pipeline has evolved from individual scripts to consolidated workflows
- Brian's original loci packets script has been forked and modified by Anh
- Additional layers of analysis have been added (Update scripts for enriching data)
- Output paths have changed to Anh's directories

## Longleaf Cluster Environment

### RHEL 9 Transition Issues (Summer 2025)
The cluster transitioned from RHEL 8 to RHEL 9, breaking several things:
- **Module versions changed**: R 4.3.1 → 4.4.0, Python 3.9.6 → 3.12.x
- **Snakemake upgraded**: Version 8.x removed the `--cluster` flag
- **Old modules unavailable**: Cannot load R 4.3.1 or Python 3.9.x anymore

### Solution: Snakemake 7 in Virtual Environment
We installed Snakemake 7.32.4 (last version with --cluster support) in a virtual environment:
```bash
# Activate the Snakemake 7 environment
source ~/snakemake7_env/bin/activate

# This version supports the old --cluster syntax!
snakemake --cluster "sbatch --time {resources.time} --mem={resources.mem_mb}" -j 8
```

### Module Management
- **Always use `module purge`** then explicitly load what you need
- **Current available versions**: R 4.4.0, Python 3.12.4
- **For reproducibility**: Specify exact versions in scripts
```bash
module purge
module load r/4.4.0
module load python/3.12.4
```

### Storage Strategy
- **/proj/raulab/**: Main project storage (good for active work and results)
- **/work/**: High IO performance for active data processing
- **/users/**: 10TB personal storage for inactive datasets
- **Symlinks work**: Can link to original genome cache to save space

## Key Commands

### Running the Pipeline (Updated for RHEL 9)
```bash
# Submit with Snakemake 7 (supports --cluster)
sbatch scripts/from_scratch/snakemake/submit_test_chr1_v2.sh

# Or run locally
source ~/snakemake7_env/bin/activate
module load r/4.4.0
snakemake --snakefile scripts/from_scratch/snakemake/Snakefile -j 2
```

## Architecture and Workflow

### Pipeline Structure (Numbered Scripts)
```
scripts/from_scratch/
├── snakemake/
│   ├── 01_ropScan.R         # Genome-wide QTL scan
│   ├── 02_permTest.R         # Permutation thresholds (20 perms for testing)
│   └── Snakefile             # Main workflow
├── 03_detectSigLoci.R        # Find significant regions
├── 04_joinLoci.R             # Merge trait & expression QTLs (miQTL only now)
├── 10_joinPositionID.R       # Consolidate unique positions
├── 11_mergeLociDetectGenes.R # Find genes in loci
├── 12_splitGenesToSpecies.R  # Split mouse/human genes
├── 20_getPathogenicity.R     # Disease associations from Open Targets
└── 21_getMouseGenePhenotypes.R # Mouse phenotypes from MouseMine
```

### Key Data Paths
- **Phenotypes**: `data/processed/phenotypes/boxcox_individual_{Ctrl,Iso}.csv`
- **Genome cache**: Symlinked to `/proj/raulab/projects/cc_gwas/data/raw/genomes/Corrected_Haplotypes_Cache_CC_83024`
- **Bulk expression**: Symlinked from `/proj/raulab/projects/cc_gwas/data/processed/joinLoci/bulk_exp/`
- **NRVM data**: Symlinked from `/proj/raulab/projects/cc_gwas/data/processed/joinLoci/nrvms/`
- **Test outputs**: `data/processed/test_chr1/` (for testing)
- **Original pipeline**: `/proj/raulab/projects/cc_gwas/` (for reference)
- **Final outputs**: `results/QTLpackets/` (loci packets with plots and Excel summaries)

### Snakemake Configuration
- **Version issue**: Must use Snakemake 7.x for --cluster support
- **Test configuration**: 2 traits (BW.day.21, LV.Mass.Corrected.21), both treatments
- **Permutations**: Reduced to 20 for testing (normally 50-100)
- **Memory allocation**: 4-6GB per job for testing

### Key R Dependencies
- **miqtl package**: From GitHub (gkeele/miqtl)
- **Core packages**: tidyverse, data.table, biomaRt
- **Analysis packages**: sva, edgeR, zoo
- **API packages**: InterMineR (MouseMine), httr (Open Targets)

## SLURM Environment
- **Job submission**: Use sbatch with resource specifications
- **Memory format**: Use `--mem=XXg` or `--mem=XXXX` (in MB)
- **Output logs**: `.slurmlogs/` directory
- **No official workflow manager support**: Run Snakemake within SLURM allocation

## Brian's Coding Style Guidelines

### General R Style
- **Libraries**: Load all libraries at the top using explicit `library()` calls
- **Command-line interfaces**: Use `optparse` for scripts requiring command-line arguments with descriptive help text
- **Piping**: Prefer tidyverse pipe `|>` over magrittr `%>%`

- **Comments**: Use section headers with `# --- Section Name ---` format for major code blocks
- **Debugging**: Include `print()` statements for parsed arguments and key progress indicators

### Data Processing Patterns
- **File I/O**: 
  - Always check/create output directories before writing: `if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }`
  - Use RDS format for R objects that need to be reloaded
  - CSV for human-readable outputs and cross-language compatibility
- **Nested lists**: Use nested named lists for organizing multi-dimensional data (e.g., `data[[trait]][[drug]]`)
- **Filtering**: Chain dplyr operations for clarity, use `filter()` with logical conditions
- **String manipulation**: Use `str_match()` with regex for pattern extraction from filenames

### Bioinformatics-Specific Patterns
- **Batch processing**: Loop through traits and treatment conditions systematically
- **Arguments and functions**: Do not specify non-standard parameters for analysis-based functions without asking for permission. In such cases, use `man` or `?FUNCTION` for bash and R commands respectively so that you can reference the documentation

### Documentation Style
- **Script headers**: Brief purpose statement at top of standalone scripts
- **Markdown notes**: Structured format with:
  - *Purpose*: What the script does
  - *Expected inputs*: Path templates with `{parameter}` notation
  - *Generated outputs*: Description of results
  - *Action items*: Known improvements or checks needed
- **Variable naming**: Descriptive snake_case for data objects, dot notation acceptable for function outputs (e.g., `miqtl.rop.scan.scaled`)

## Environment Setup (As of Aug 2025)

### Consolidated miqtl-env
All dependencies are now in a single conda environment:
- **Location**: `~/mambaforge/envs/miqtl-env/`
- **Key packages**: R 4.4.0, Python 3.11, Snakemake 7.32.4, tidyverse, biomaRt, DESeq2, edgeR, sva
- **Environment spec**: `envs/environment.yml`
- **Rebuild script**: `scripts/rebuild_env.sh`
- **miqtl package**: Installed from GitHub to `/proj/raulab/users/brian/cc_gwas/R_packages/`

### Submission & Logging
- **Main script**: `scripts/from_scratch/snakemake/submit_with_existing.sh`
- **Log structure**: `.slurmlogs/YYYY-MM-DD/HHMMSS_JOBID/` with organized subdirectories
- **Job naming**: Cluster jobs named by Snakemake rule for clarity

## Known Issues & Solutions

### Snakemake --cluster Flag Error
**Problem**: Snakemake 8.x removed --cluster support
**Solution**: miqtl-env has Snakemake 7.32.4 installed with mamba

### Missing R/Python Versions
**Problem**: R 4.3.1 and Python 3.9.6 no longer available after RHEL 9 transition
**Solution**: Using R 4.4.0 and Python 3.11 in miqtl-env

### InterMineR/RCurl Library Issues
**Problem**: InterMineR fails with missing ICU libraries
**Status**: InterMineR installed from GitHub but has shared library dependencies
**Current Workaround**: Hardcoded LD_LIBRARY_PATH in submit script to /nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib
**TODO**: Find a solution that doesn't involve hardcoding the library path - should dynamically detect conda env lib path

### pyLMM Results
**Status**: COMPLETELY REMOVED from pipeline - focusing on miQTL only
**Scripts affected**: 
- 04_joinLoci.R - Modified to skip pyLMM processing
- 13_multTrait_cis-eQTL_nrvmExp.R - Created modified version without PyLMM cis-eQTL dependencies
- 22_makeLociPackets.R - Created modified version with PyLMM plotting sections commented out

### MouseMine Database Access
**Problem**: MouseMine appears to be offline/retired as of Aug 29, 2025
**Impact**: Script 21_getMouseGenePhenotypes.R cannot retrieve phenotype annotations
**Available mines**: AllianceMine, AquaMine, FlyMine, HumanMine, HymenopteraMine, MaizeMine, PhytoMine, PlanMine, TargetMine, ThaleMine, WormMine
**TODO**: Find alternative source for mouse phenotype annotations or skip this enrichment step
- Avoid making new files to test things. Instead, run bash commands directly or change the file that you're testing (while using git for easy managemetn of versions)