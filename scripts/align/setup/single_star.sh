#!/bin/bash

#SBATCH --job-name=star_align_CC_plate_1
#SBATCH --output=.slurmlogs/%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=380G
#SBATCH --time=48:00:00

module load star/2.7.9a

sample="CC_Plate_1"
r1="data/processed/merged_fastq/${sample}_R1.fastq.gz"
r2="data/processed/merged_fastq/${sample}_R2.fastq.gz"
barcodes=data/raw/barcodes/v5B_96_clean.txt
outdir="data/processed/aligned/${sample}"

STAR --runMode alignReads \
    --outSAMmapqUnique 60 \
    --runThreadN 24 \
    --outSAMunmapped Within \
    --soloStrand Forward \
    --quantMode GeneCounts \
    --outBAMsortingThreadN 8 \
    --genomeDir data/processed/star_index/ \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 \
    --soloCBlen 14 \
    --soloUMIstart 15 \
    --soloUMIlen 14 \
    --soloCellFilter None \
    --soloUMIdedup NoDedup \
    --soloCBwhitelist ${barcodes} \
    --soloFeatures Gene \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outFilterMultimapNmax 1 \
    --readFilesCommand zcat \
    --limitBAMsortRAM 35004580355 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${outdir}/ \
    --readFilesIn ${r2} ${r1} \
    --outTmpDir /work/users/b/g/bgural #This needs to be changed if someone else runs it. 