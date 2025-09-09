#!/bin/bash

# Script to set up Snakemake 7.x in a virtual environment
# This version still supports the --cluster flag

echo "Setting up Snakemake 7.x environment..."

# Load Python module
module purge
module load python/3.12.4

# Create virtual environment
python -m venv ~/snakemake7_env

# Activate it
source ~/snakemake7_env/bin/activate

# Install Snakemake 7.32.4 (last 7.x version)
pip install --upgrade pip
pip install snakemake==7.32.4
pip install pandas pyyaml

echo ""
echo "Snakemake 7.x environment created!"
echo "Installed version:"
snakemake --version
echo ""
echo "To use this environment, add this to your scripts:"
echo "  source ~/snakemake7_env/bin/activate"
echo ""
echo "This version supports the --cluster flag!"