#!/bin/bash

#SBATCH --job-name=snake1
#SBATCH --time=2:00:00
#SBATCH --mem=4000
#SBATCH --output=.slurmlogs/slurmlog1/%j.out


module load r
module load python

snakemake --snakefile Scripts/Snakefile_part1 --forcerun all -j 1000 --latency-wait 60 --cluster "sbatch --time 1:00:00 --mem=2000 -N 1 -n 1 -o ./.slurmlogs/slurm1/%j.out" --keep-going


snakemake --rulegraph | dot -Tsvg > repo/dag.svg

#run this line in terminal: snakemake -s Scripts/Snakefile_part1 --unlock
#run this line in terminal: sbatch Scripts/snake_1.sh
#to check how many files are done: ls -1 someFolderPath | wc -l