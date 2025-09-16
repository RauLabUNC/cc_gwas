# --- Load Libraries ---
library(miqtl)
library(tidyverse)
library(optparse)
source("scripts/helpers/scan_h2lmm.R")

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
# Genome Scan Setup
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Read the processed phenotype file
phenotypes <- read.csv(opt$input, check.names = FALSE)
colnames(phenotypes)
# Ensure covariate type
if (!"Sex" %in% names(phenotypes)) stop("Sex column missing in phenotypes")
phenotypes$Sex <- as.factor(phenotypes$Sex)

qtl_trait <- opt$qtl_trait
print(qtl_trait)

# Verify trait column and build robust formula (handles special chars via backticks)
if (!qtl_trait %in% names(phenotypes)) {
  stop(sprintf("Trait '%s' not found in input data", qtl_trait))
}
form <- stats::as.formula(paste0("`", qtl_trait, "` ~ 0 + Sex"))

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
  formula = form,  # do not use get(qtl_trait)
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