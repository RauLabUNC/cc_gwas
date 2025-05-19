#!/bin/bash
#SBATCH --job-name=unzip_haplo
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --output=unzip_haplo_%j.out
#SBATCH --error=unzip_haplo_%j.err

# Load any necessary modules
module load unzip

# Define file paths
ZIP_FILE="/proj/raulab/projects/cc_gwas/data/raw/genomes/Corrected_Haplotypes_Cache_CC_83024.zip"
TARGET_DIR="/proj/raulab/projects/cc_gwas/data/raw/genomes"

# Navigate to the target directory
cd ${TARGET_DIR}

# Remove the existing directory
rm -rf haplotype_cache_cc_083024

# Unzip the file
unzip "${ZIP_FILE}"

echo "Unzipping completed successfully" 