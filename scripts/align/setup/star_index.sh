#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=24:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=8
#SBATCH --output=.slurmlogs/%j.out

# Load STAR module
module load star/2.7.9a

# Set directories 
GENOME_DIR="data/raw/reference"
INDEX_DIR="data/processed/star_index"

# Create index directory
mkdir -p $INDEX_DIR

# Generate STAR index
STAR --runMode genomeGenerate \
    --genomeDir $INDEX_DIR \
    --genomeFastaFiles $GENOME_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --sjdbGTFfile $GENOME_DIR/Mus_musculus.GRCm39.110.gtf \
    --runThreadN 8