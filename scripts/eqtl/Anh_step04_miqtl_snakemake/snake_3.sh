#!/bin/bash

#SBATCH --job-name=snake3
#SBATCH --time=01:00:00
#SBATCH --mem=1000
#SBATCH --output=.slurmlogs/slurmlog3/%j.out


module load r
module load python

#snakemake --snakefile Scripts/Snakefile_part3 --forcerun all -j 1000 --latency-wait 60 --cluster "sbatch --time 00:01:00 --mem=100 -N 1 -n 1 -o ./.slurmlogs/slurm3/%j.out" --keep-going

snakemake --snakefile Scripts/Snakefile_part3_allSNPs_1Chr --forcerun all -j 1000 --latency-wait 60 --cluster "sbatch --time 00:01:00 --mem=100 -N 1 -n 1 -o ./.slurmlogs3/%j.out" --keep-going


snakemake --rulegraph | dot -Tsvg > repo/dag.svg

#run this line in terminal: snakemake -s Scripts/Snakefile_part3 --unlock
#or     snakemake -s Scripts/Snakefile_part3_allSNPs_1Chr --unlock
#run this line in terminal: sbatch Scripts/snake_3.sh