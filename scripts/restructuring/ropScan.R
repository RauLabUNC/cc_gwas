# --- Load Libraries ---
library(miqtl)
library(tidyverse)
library(optparse)
source("scripts/extras/scan_h2lmm.R")

# --- Define Command-Line Arguments ---
option_list <- list(
  make_option(c("--input"), type = "character", help = "Path to input processed phenotype file"),
  make_option(c("--output"), type = "character", help = "Path to output RDS file"),
  make_option(c("--normalization"), type = "character", help = "Normalization method (e.g., zscore, boxcox)"),
  make_option(c("--aggregation"), type = "character", help = "Aggregation method (e.g., individual, mean)"),
  make_option(c("--qtl_trait"), type = "character", help = "QTL trait to analyze"),
  make_option(c("--drug"), type = "character", help = "Drug treatment (e.g., Ctrl, Iso)")
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

# Run genome scan
miqtl.rop.scan.scaled <- scan.h2lmm.test(
  genomecache = genomecache,
  data = phenotypes,
  pheno.id = "gwas_temp_id",
  geno.id = "Strain",
  formula = get(qtl_trait) ~ 0 + Sex,  # Incorporate covariates
  use.multi.impute = FALSE,
  return.allele.effects = TRUE,
  use.fix.par = TRUE
)

# --- Save Output ---
# Ensure the output directory exists
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the output as an RDS file
saveRDS(miqtl.rop.scan.scaled, opt$output)