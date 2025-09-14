#!/bin/bash
#SBATCH --job-name=merge_p${1}
#SBATCH --time=1:00:00
#SBATCH --mem=4000
#SBATCH --output=.slurmlogs/%j.out

# Create output directory if it doesn't exist
mkdir -p data/processed/merged_fastq

# Select appropriate lanes based on plate number
if [[ ${1} == "3" || ${1} == "4" ]]; then
    # Plates 3 and 4 only use lane 5
    cat data/rnaseq/CC_Plate_${1}_S${2}_L005_R1_001.fastq.gz > data/processed/merged_fastq/CC_Plate_${1}_R1.fastq.gz
    cat data/rnaseq/CC_Plate_${1}_S${2}_L005_R2_001.fastq.gz > data/processed/merged_fastq/CC_Plate_${1}_R2.fastq.gz
elif [[ ${1} == "5" ]]; then
    # Plate 5 uses lanes 1-4
    cat data/rnaseq/CC_Plate_${1}_S${2}_L00{1,2,3,4}_R1_001.fastq.gz > data/processed/merged_fastq/CC_Plate_${1}_R1.fastq.gz
    cat data/rnaseq/CC_Plate_${1}_S${2}_L00{1,2,3,4}_R2_001.fastq.gz > data/processed/merged_fastq/CC_Plate_${1}_R2.fastq.gz
else
    # Plates 1 and 2 use lanes 1-4,6-8
    cat data/rnaseq/CC_Plate_${1}_S${2}_L00{1,2,3,4,6,7,8}_R1_001.fastq.gz > data/processed/merged_fastq/CC_Plate_${1}_R1.fastq.gz
    cat data/rnaseq/CC_Plate_${1}_S${2}_L00{1,2,3,4,6,7,8}_R2_001.fastq.gz > data/processed/merged_fastq/CC_Plate_${1}_R2.fastq.gz
fi