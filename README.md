# Collaborative Cross Heart Failure GWAS Pipeline

This repository contains the computational pipeline for identifying genetic and transcriptomic determinants of heart failure susceptibility in Collaborative Cross mice, as described in Kimball et al. (2026) DOI: ####.

## Overview

We induced isoproterenol-mediated heart failure in 63 diverse Collaborative Cross strains (~400 mice total) to map genetic variants controlling cardiac remodeling response. This pipeline integrates:
- Trait QTL mapping - identifying genomic loci associated with cardiac phenotypes
- RNA-seq alignment - processing bulk RNA-seq from cardiac tissue  
- eQTL mapping - linking genetic variation to gene expression
- Gene prioritization - automated evidence packets for functional validation

## Pipeline Structure

The workflow consists of four semi-independent components with the following dependencies:

```
align (RNA-seq) ────────┐
                        ├──> eQTL ──┐
traitqtl (QTL) ────────┤            ├──> packets (gene prioritization)
                       └────────────┘
```

### 1. Trait QTL Mapping (`scripts/traitqtl/`)
**Status:** Fully integrated  
**Dependencies:** None (requires only phenotype data and genomecache)

Maps genetic loci associated with cardiac phenotypes using multi-parent regression-on-probabilities. Performs Box-Cox normalization, genome scans, and permutation-based significance testing.

- **Entry point:** `sbatch scripts/traitqtl/snakemake/snake.sh`
- **Outputs:** 
  - `data/processed/joinLoci/relational_tables/genesInLoci.rds` - genes within QTL intervals
  - Significance thresholds and scan results

### 2. RNA-seq Alignment (`scripts/align/`)
**Status:** Works independently  
**Dependencies:** None (requires raw FASTQ files)

STAR-based alignment pipeline for BRB-seq libraries with 14bp sample barcodes and UMIs. Generates count matrices used for eQTL mapping and differential expression.

- **Entry point:** See `scripts/align/README.md` for execution details
- **Outputs:**
  - `data/processed/star_aligned/CC_Plate_*/` - aligned BAMs and count matrices
  - `data/processed/expression/cc_raw_counts.csv` - merged count matrix

### 3. eQTL Mapping (`scripts/eqtl/`)
**Status:** Partially integrated, needs work  
**Dependencies:** 
  - Required: RNA-seq alignment outputs
  - Optional: trait QTL results (to restrict analysis to genes in trait loci)

Maps expression QTLs for genes, with optional restriction to those within trait QTL intervals. Currently split across three separate snake scripts with hardcoded paths requiring manual execution.

- **Entry point:** Currently broken - see `scripts/eqtl/README.md` for details
- **Integration status:** Needs consolidation into single snake.sh with conda environment

### 4. Gene Packet Generation (`scripts/packets/`)
**Status:** Core pipeline works, needs manual post-processing  
**Dependencies:** All three upstream components

Combines trait QTL, eQTL, and pathogenicity data into standardized "gene packets" - evidence summaries for each locus including:
- LocusZoom-style visualizations
- Founder mutation tables
- Multi-sheet Excel workbooks with expression, disease associations, and phenotypes
- Automated cardiac-relevance flagging

- **Entry point:** `sbatch scripts/packets/snakemake/snake.sh`
- **Outputs:** `results/qtl_packets/` - one packet per locus cluster

**Note:** Manual curation scripts (`scripts/packets/Anh_*.R`) for adding CPM values, eQTL tags, and protein annotations need integration into Snakemake workflow.

## Quick Start

### Prerequisites

```bash
# Activate consolidated conda environment
set +u
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env
set -u
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
```

Environment contains:
- R 4.4.0 with miQTL, DESeq2, biomaRt, and related packages
- Python 3.11
- Snakemake 7.32.4

### Required Input Data

1. **Genotypes:** Collaborative Cross founder haplotypes
   - Location: `data/raw/genomes/haplotype_cache_cc_083024/`
   - Source: [UNC Systems Genetics](https://csbio.unc.edu/CCstatus/CCGenomes/) for raw files, then conver to HAPPY format

2. **Phenotypes:** Processed cardiac measurements
   - Will be deposited to data repository (upload in progress)
   - Traits: body weight, chamber dimensions, wall thickness, fibrosis, cross-sectional area

3. **RNA-seq:** BRB-seq FASTQ files (paired-end)
   - 5 libraries, ~440 samples with 14bp barcodes
   - Available in SRA (see Data Availability)

4. **Reference genome:** GRCm39 (Ensembl 110)
   - FASTA and GTF files in `data/raw/reference/`

### Running the Pipeline

**Test mode** (Chr 1 only, 2 traits, 20 permutations):
```bash
export SNAKEMAKE_MODE="test"
sbatch scripts/traitqtl/snakemake/snake.sh
```

**Full mode** (all chromosomes, 10 traits, 100 permutations):
```bash
export SNAKEMAKE_MODE="full"  # or omit, as "full" is default
sbatch scripts/traitqtl/snakemake/snake.sh
```

For complete workflow execution order, see component-specific READMEs.

## Pipeline Configuration

All Snakemake workflows use `config.yaml` files in their respective directories. Test vs. full mode is controlled by the `SNAKEMAKE_MODE` environment variable, which determines:
- Number of chromosomes scanned (1 vs all)
- Traits analyzed (2 vs 10)
- Permutation count (20 vs 100)
- Output file prefixes (`test_` prefix in test mode)

Modify resource allocations (time, memory) in the `resources:` section of each `config.yaml`.

## Logging

All Snakemake workflows use standardized date-based logging:
- Format: `.slurmlogs/YYYY-MM-DD/HHMMSS_JOBID/`
- Contains:
  - `main.out/err` - controller job logs
  - `jobs/RULE/WILDCARDS/*.out/err` - individual rule logs
  - `snakemake_stats.json` - job statistics

## Current Limitations
This pipeline, although documented, lacks complete integration and streamlining.
1. **eQTL pipeline fragmentation** - Split across 3 snake scripts requiring sequential manual execution. Needs consolidation into single entry point.

2. **Hardcoded paths** - Some scripts contain absolute paths to `/proj/raulab/users/Anh/`. Need parameterization for portability.

3. **Manual post-processing** - `scripts/packets/Anh_*.R` scripts for CPM addition, eQTL tagging, and protein annotations are not integrated into Snakemake.

4. **Incomplete integration** - Four components can run independently but lack unified master workflow.

## Repository Structure

```
.
├── data/
│   ├── raw/
│   │   ├── genomes/           # CC haplotype cache
│   │   ├── phenotypes/        # Raw cardiac measurements  
│   │   ├── rnaseq/            # FASTQ files
│   │   ├── reference/         # GRCm39 genome + GTF
│   │   └── barcodes/          # BRB-seq barcodes
│   └── processed/             # Intermediate files
├── results/                   # Final outputs
│   └── qtl_packets/          # Gene evidence packets
├── scripts/
│   ├── align/                # RNA-seq pipeline
│   ├── traitqtl/             # QTL mapping
│   ├── eqtl/                 # eQTL mapping  
│   ├── packets/              # Packet generation
│   ├── figures/              # Standalone plotting scripts
│   └── helpers/              # Utility functions
├── envs/
│   └── environment.yml       # Conda environment spec
└── README.md
```

## Data Availability

**RNA-seq data:**
- NRVM data: Sequence Read Archive, BioProject ID 1328793
- CC mice cardiac tissue: Sequence Read Archive, BioProject ID 1331689

**Phenotype data:** To be deposited in data repository (upload scheduled)



## Citation

Kimball et al. (2026). "Genetic and transcriptomic determinants of heart failure susceptibility in Collaborative Cross mice." Journal TBD. DOI: ####

For the gene prioritization methodology, see companion publication:
Gural et al. (in preparation). "Systematic evidence integration for gene prioritization from large QTL intervals in Collaborative Cross mice." BMC Research Notes.

Rau Lab, University of North Carolina at Chapel Hill
