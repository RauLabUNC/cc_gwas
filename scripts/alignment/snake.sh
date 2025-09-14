#!/bin/bash

#SBATCH --job-name=snake
#SBATCH --time=120:00:00
#SBATCH --mem=4G
#SBATCH --output=.slurmlogs/%j.out

module load python/3.9.6
module load star/2.7.9
module load fastqc/0.12.1 
module load seqkit/2.8.2
module load salmon
module load cutadapt

snakemake --forcerun all -j 50 --latency-wait 60 \
    --cluster "sbatch \
        --time=120:00:00 \
        --mem={resources.mem_mb}M \
        -c {threads} \
        -N 1 \
        -o ./.slurmlogs/%j.out" \
    --keep-going