# Load libraries
library(miqtl)
library(tidyverse)
# Read in the first trailing argument for phenotype and treatment
args <- commandArgs(trailingOnly = TRUE)

if(args[2] == "control"){
  treatment <- "0"
}else{
  treatment <- "1"
}

# Load the data
scan_file <- file.path("data/processed/scans", args[2], paste0(args[1], "_scan_results.rds"))
scan <- readRDS(scan_file)

# Load phenotypes
phenotypes <- read.csv("data/processed/phenotypes/meanCenterScaledByTreat_03242025.csv")

# Remove NAs, summarize groups to means, and scale
scaled_phenotypes <- phenotypes |> 
  filter(!is.na(args[1]) & Drug_Binary == treatment)

# Load genome cache
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Permute phenotypes and run scans on each
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan, 
                                                      method = "permutation", num.samples = 50,
                                                      use.BLUP = T, model.type = "null")

permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = scaled_phenotypes,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1)

permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, 
                                        percentile = 0.85) # 10.1186/s40168-023-01552-8 uses 85 percentile


# Ensure the directory exists
output_dir <- file.path("data/processed/scan_thresholds", args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_threshold <- file.path(output_dir, paste0(args[[1]], "_threshold.rds"))
output_scan <- file.path(output_dir, paste0(args[[1]], "_scan.rds"))

# Save the threshold and scan as RDS files
saveRDS(permute_threshold, output_threshold)
saveRDS(permuted_scans, output_scan)
