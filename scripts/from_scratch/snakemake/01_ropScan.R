# --- Load Libraries ---
library(miqtl)
library(tidyverse)
library(optparse)
source("scripts/extras/scan_h2lmm.R")

# --- Define Command-Line Arguments ---
option_list <- list(
  make_option(c("--input"), type = "character", help = "Path to input processed phenotype file"),
  make_option(c("--output_scan"), type = "character", help = "Path to output RDS file of scan"),
  make_option(c("--output_perms"), type = "character", help = "Path to output RDS file of permuted phenotypes"),
  make_option(c("--normalization"), type = "character", help = "Normalization method (e.g., zscore, boxcox)"),
  make_option(c("--aggregation"), type = "character", help = "Aggregation method (e.g., individual, mean)"),
  make_option(c("--qtl_trait"), type = "character", help = "QTL trait to analyze"),
  make_option(c("--drug"), type = "character", help = "Drug treatment (e.g., Ctrl, Iso)"),
  make_option(c("--mode"), type = "character", default = "full", help = "Run mode: 'full' or 'test' (test runs chr1 only)")
)

# Parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Print the parsed arguments (for debugging)
print(opt)

# --- Load Input Data ---
# Read the processed phenotype file
phenotypes <- read.csv(opt$input)

# --- Genome Scan Setup ---
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

qtl_trait <- opt$qtl_trait
print(qtl_trait)

# Set chromosome parameter based on mode
if (opt$mode == "test") {
  cat("Running in TEST mode - chromosome 1 only\n")
  chr_to_scan <- 1
} else {
  cat("Running in FULL mode - all chromosomes\n")
  chr_to_scan <- "all"
}

# Run genome scan
miqtl.rop.scan.scaled <- scan.h2lmm.test(
  genomecache = genomecache,
  data = phenotypes,
  pheno.id = "gwas_temp_id",
  geno.id = "Strain",
  formula = get(qtl_trait) ~ 0 + Sex,  # Incorporate covariates
  use.multi.impute = FALSE,
  return.allele.effects = TRUE,
  use.fix.par = TRUE,
  chr = chr_to_scan
)

# Generate permuted phenotypes from the scan object

num_perms <- ifelse(opt$mode == "test", 50, 1000)
cat(sprintf("\nGenerating %d permutations...\n", num_perms))
permuted_phenotype <- generate.sample.outcomes.matrix(
  scan.object = miqtl.rop.scan.scaled,
  method = "permutation",
  num.samples = num_perms,
  use.BLUP = TRUE,
  model.type = "null"
)

# --- Save Output ---
# Ensure the output directory exists
output_scan_dir <- dirname(opt$output_scan)
if (!dir.exists(output_scan_dir)) {
  dir.create(output_scan_dir, recursive = TRUE)
}
output_perms_dir <- dirname(opt$output_perms)
if (!dir.exists(output_perms_dir)) {
  dir.create(output_perms_dir, recursive = TRUE)
}


# Save the output as an RDS file
saveRDS(miqtl.rop.scan.scaled, opt$output_scan)
saveRDS(permuted_phenotype, opt$output_perms)