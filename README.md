# Collaborative Cross GWAS/QTL Mapping Project

This repository contains the genetics analysis code and data for the collaborative cross heart failure project conducted by the Rau Laboratory at the University of North Carolina at Chapel Hill (UNC-CH).

## Introduction

The Collaborative Cross (CC) is a powerful genetic resource that enables the study of complex traits in mice. It is a panel of recombinant inbred strains derived from eight founder strains, designed to capture the genetic diversity found in natural populations. In this project, we leverage the CC to perform Genome-Wide Association Studies (GWAS) and Quantitative Trait Locus (QTL) mapping to identify genetic variants associated with heart failure. 

![CC structure](https://github.com/RauLabUNC/cc_gwas/blob/main/repo/gr1_lrg.jpg?raw=true)

**Collaborative Cross: A tool for discovery population genetics at scale**

This figure illustrates the genetic diversity and structure of the mouse strains used in the study. The founder strains consist of five common laboratory inbred strains and three wild-derived inbred strains, capturing a significant portion of the genetic variation observed in Mus musculus. The Collaborative Cross strains are eight-way recombinant inbred strains derived from the founders, while the Diversity Outbred mice (not used in this study) are the result of continuous outbreeding, leading to dense recombination of genome structure and high heterozygosity. These mouse strains serve as valuable tools for investigating the genetic basis of complex traits. Figure sourced from [Saul 2019](https://www.cell.com/trends/genetics/fulltext/S0168-9525%2819%2930065-4).

## Repository Structure

- `data/`: Contains the raw and processed data used in the analysis.
- `scripts/`: Contains the R scripts used for data preprocessing, statistical analysis, and visualization.
- `results/`: Contains the output files and figures generated from the analysis.

## Getting Started

Once the repo is completed, users may reproduce the analysis by doing the following:

1. Clone this repository to a high-performance computing cluster with SLURM job scheduling.
2. Install the required dependencies listed in the `requirements.txt` file.*
3. Run `bash snake.sh` in the main project directory to execute the Snakemake job scheduler

## Snakemake Workflow

The order of script execution is visualized in this directed acyclic graph:

![DAG Plot](https://github.com/RauLabUNC/cc_gwas/blob/main/repo/dag.svg?raw=true)

## To-Do
This repo is still in its early stages, with many remaining items to complete.
- [ ] Integrate GWAS and QTL mapping pipelines for physiological traits
- [ ] Investigate the functional implications of the identified genetic variants.
- [ ] Write manuscript summarizing the findings and submit it for publication.

For any questions or inquiries, please contact Christoph Rau.
