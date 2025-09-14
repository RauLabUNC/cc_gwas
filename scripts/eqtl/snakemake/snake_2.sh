#!/bin/bash

#SBATCH --job-name=snake2
#SBATCH --time=2:00:00
#SBATCH --mem=4000
#SBATCH --output=.slurmlogs/slurmlog2/%j.out


module load r
module load python

snakemake --snakefile Scripts/Snakefile_part2_oneChr --forcerun all -j 1000 --latency-wait 60 --cluster "sbatch --time 1:15:00 --mem=1200 -N 1 -n 1 -o ./.slurmlogs/slurm2/%j.out" --keep-going

snakemake --rulegraph | dot -Tsvg > repo/dag.svg


#run this line in terminal: snakemake -s Scripts/Snakefile_part2_oneChr --unlock
#run this line in terminal: sbatch Scripts/snake_2.sh

#to check how many files are done: ls -1 /proj/raulab/users/Anh/CC-miQTL_VST/Data/Scan_thresholds/Permute_1Gene1Chr/250602-VST_614/iso | wc -l