# This script is meant to plot genome scans by miQTL with their permuted...
# ...significance thresholds

# Load libraries
library(miqtl)

# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]


## Load the data
# Scan
scan_file <- file.path("data/processed/scans", args[2], paste0(as.character(phenotype_of_interest), "_scan_results.rds"))
scan <- readRDS(scan_file)

# Threshold
threshold_file <- file.path("data/processed/scan_thresholds", args[2], paste0(as.character(phenotype_of_interest), "_threshold.rds"))
threshold <- readRDS(threshold_file)

# Plot the results
output_dir <- "results/genome_scans_thresholds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

png_name <- file.path(output_dir, args[2], paste0(as.character(phenotype_of_interest), "_scan_threshold.png"))
png(file = png_name,
    width = 8, 
    height = 4,
    units = "in",
    res = 600)

genome.plotter.whole(scan.list=list(ROP = scan), 
                     hard.thresholds = threshold)

dev.off()
