#!/usr/bin/env Rscript
# Minimal test script for QTL scanning with a single gene/trait
# This tests the basic pipeline with correct paths

cat("====================================\n")
cat("Testing Single QTL Scan\n")
cat("====================================\n\n")

# Load libraries
library(miqtl)
library(tidyverse)

# Source the helper script (from original location)
source("/proj/raulab/projects/cc_gwas/scripts/extras/scan_h2lmm.R")

# Define paths - pointing to ORIGINAL data locations
genomecache <- "/proj/raulab/projects/cc_gwas/data/raw/genomes/Corrected_Haplotypes_Cache_CC_83024"
phenotype_file <- "/proj/raulab/projects/cc_gwas/data/processed/phenotypes/meanCenterScaledByTreat_03242025.csv"

# Check if files exist
cat("Checking file paths:\n")
cat("  Genome cache: ", file.exists(genomecache), "\n")
cat("  Phenotype file: ", file.exists(phenotype_file), "\n")
cat("  Helper script: ", file.exists("/proj/raulab/projects/cc_gwas/scripts/extras/scan_h2lmm.R"), "\n\n")

# Load phenotype data
cat("Loading phenotype data...\n")
phenotypes <- read.csv(phenotype_file)

# Check structure
cat("Phenotype data dimensions: ", nrow(phenotypes), "x", ncol(phenotypes), "\n")
cat("Column names (first 10):\n")
print(head(names(phenotypes), 10))
cat("\n")

# Check for required columns and fix naming
required_cols <- c("gwas_temp_id", "Strain", "Sex")
if(!"Sex" %in% names(phenotypes) && "Sex_Binary" %in% names(phenotypes)) {
  phenotypes$Sex <- phenotypes$Sex_Binary
  cat("NOTE: Renamed Sex_Binary to Sex\n")
}
missing_cols <- setdiff(required_cols, names(phenotypes))
if(length(missing_cols) > 0) {
  cat("WARNING: Missing required columns: ", paste(missing_cols, collapse=", "), "\n")
  cat("Available columns with 'sex':\n")
  print(names(phenotypes)[grep("sex|Sex", names(phenotypes), ignore.case=TRUE)])
}

# Pick a test trait - use BW.day.0 which should be a real phenotype
numeric_cols <- names(phenotypes)[sapply(phenotypes, is.numeric)]
# Remove obvious ID columns and binary columns
trait_cols <- numeric_cols[!grepl("id|ID|index|Binary|Sex", numeric_cols, ignore.case=TRUE)]
# Use BW.day.0 if available, otherwise first available trait
test_trait <- if("BW.day.0" %in% trait_cols) "BW.day.0" else trait_cols[1]

cat("\nUsing test trait: ", test_trait, "\n")
cat("Trait summary:\n")
print(summary(phenotypes[[test_trait]]))

# Test with small chromosome first
cat("\n====================================\n")
cat("Running QTL scan on chromosome 19 (smallest)\n")
cat("====================================\n\n")

# Run scan with minimal settings
tryCatch({
  miqtl_scan <- scan.h2lmm.test(
    genomecache = genomecache,
    data = phenotypes,
    pheno.id = "gwas_temp_id",  # Adjust if different column name
    geno.id = "Strain",
    formula = as.formula(paste(test_trait, "~ 0 + Sex")),
    use.multi.impute = FALSE,
    return.allele.effects = TRUE,
    use.fix.par = TRUE,
    chr = "19"  # Just chromosome 19 for testing
  )
  
  cat("SUCCESS: Scan completed!\n")
  cat("Result structure:\n")
  print(str(miqtl_scan, max.level = 2))
  
  # Save output
  output_file <- paste0("test_scan_", test_trait, "_chr19.rds")
  saveRDS(miqtl_scan, output_file)
  cat("\nResults saved to: ", output_file, "\n")
  
}, error = function(e) {
  cat("ERROR during scan:\n")
  cat(e$message, "\n\n")
  cat("Debugging info:\n")
  cat("  Unique strains: ", length(unique(phenotypes$Strain)), "\n")
  cat("  Sex values: ", unique(phenotypes$Sex), "\n")
  cat("  Missing values in trait: ", sum(is.na(phenotypes[[test_trait]])), "\n")
})