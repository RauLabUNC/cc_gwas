#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --time=120:00:00
#SBATCH --mem=1000
#SBATCH --output=.slurmlogs/%j.out

module load r r/4.3.1
module load python/3.9.6

snakemake -j 50 --latency-wait 60 --cluster "sbatch --time 120:00:00 --mem=3000 -N 1 -n 1 -o ./.slurmlogs/%j.out" --keep-going

snakemake --rulegraph | dot -Tsvg > repo/dag.svg