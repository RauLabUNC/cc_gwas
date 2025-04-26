# This script is meant to follow up the genome scans and threshold testing
# It pulls out any haplotypes above the significance threshold and merges

# Load libraries
library(miqtl)
library(tidyverse)

args <- c("CSA", "iso")
# Load arguments
args <- commandArgs(trailingOnly = TRUE)

# Load the data
scan_file <- file.path("data/processed/scans", args[2], paste0(args[1], "_scan_results.rds"))
scan <- readRDS(scan_file)

threshold_file <- file.path("data/processed/scan_thresholds", args[2], paste0(args[1], "_threshold.rds"))
threshold <- readRDS(threshold_file)

threshold <- 10^(-threshold) # convert from -log10

# pull out important info for each sig loci
sig.loci <- scan$p.value[scan$p.value < threshold]

# Save if any of the 
if(length(sig.loci) > 0){
  chromosomes <- scan$chr[scan$p.value < threshold]
  positions <- scan$pos$Mb[scan$p.value < threshold]
  
  # combine them in a single df
  sig.loci.df <- data.frame(loci = names(sig.loci), p.value = sig.loci, 
                            phenotype = args[1], treatment = args[2],
                            chromosome = chromosomes, position = positions, threshold = threshold)
  }else{ #Make empty df if no loci are significant
  sig.loci.df <- data.frame(loci = NA, p.value = NA, 
                            phenotype = args[1], treatment = args[2],
                            chromosome = NA, position = NA, threshold = threshold)
}

# Ensure the directory exists
output_dir <- file.path("data/processed/sig_loci", args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
  
# Define the output file path based on the phenotype of interest
loci.file <- file.path(output_dir, 
                       paste0(args[1], "_loci.csv"))
  
# Save the loci as csv file
write.csv(sig.loci.df, loci.file, row.names = F)

