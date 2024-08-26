#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=240
#SBATCH --mem=1000

module load r r/4.3.1
module load python/3.9.6


snakemake -j 36 --latency-wait 60 --cluster "sbatch --mem=3000 -N 1 -n 1 -o ./.slurmlogs/%j.out" --keep-going

snakemake --rulegraph | dot -Tsvg > dag.svg