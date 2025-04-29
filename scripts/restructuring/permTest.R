# --- Load Libraries ---
library(miqtl)
library(tidyverse)
library(optparse)

# --- Define Command-Line Arguments ---
option_list <- list(
  make_option(c("--input_scan"), type = "character", help = "Path to input scan results RDS file"),
  make_option(c("--input_pheno"), type = "character", help = "Path to input processed phenotype file"),
  make_option(c("--output_threshold"), type = "character", help = "Path to output threshold RDS file"),
  make_option(c("--output_scan"), type = "character", help = "Path to output permuted scan RDS file"),
  make_option(c("--qtl_trait"), type = "character", help = "QTL trait to analyze"),
  make_option(c("--drug"), type = "character", help = "Drug treatment (e.g., Ctrl, Iso)")
)

# Parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Print the parsed arguments (for debugging)
print(opt)


# --- Load Input Data ---
# Load the scan results
scan <- readRDS(opt$input_scan)

# Load the processed phenotype file
phenotypes <- read.csv(opt$input_pheno)

# Filter the data based on the drug treatment
phenotypes <- phenotypes |> filter(!is.na(opt$qtl_trait) & Drug == opt$drug)

# --- Genome Scan Setup ---
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Permute phenotypes and run scans on each
permuted_phenotype <- generate.sample.outcomes.matrix(
  scan.object = scan, 
  method = "permutation", 
  num.samples = 30,
  use.BLUP = TRUE, 
  model.type = "null"
)

permuted_scans <- run.threshold.scans(
  sim.threshold.object = permuted_phenotype, 
  keep.full.scans = TRUE,
  genomecache = genomecache, 
  data = phenotypes,
  use.multi.impute = FALSE, 
  scan.seed = 1
)

permute_threshold <- get.gev.thresholds(
  threshold.scans = permuted_scans, 
  percentile = 0.9
)

# --- Save Outputs ---
# Ensure the output directory exists
output_dir <- dirname(opt$output_threshold)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the threshold and scan as RDS files
saveRDS(permute_threshold, opt$output_threshold)
saveRDS(permuted_scans, opt$output_scan)