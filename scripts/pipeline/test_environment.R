#!/usr/bin/env Rscript
# Test script to verify all required packages are available

cat("====================================\n")
cat("Testing CC-miQTL R Environment\n")
cat("====================================\n\n")

# Function to test package loading
test_package <- function(pkg_name, github_repo = NULL) {
  cat(sprintf("Testing %-30s ... ", pkg_name))
  
  tryCatch({
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
    cat("✓ LOADED\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ FAILED\n")
    if (!is.null(github_repo)) {
      cat(sprintf("  Install with: remotes::install_github('%s')\n", github_repo))
    } else {
      cat(sprintf("  Install with: install.packages('%s')\n", pkg_name))
    }
    return(FALSE)
  })
}

# Core packages
cat("Core R Packages:\n")
cat("----------------\n")
core_packages <- c("tidyverse", "dplyr", "tidyr", "ggplot2", "readr", 
                  "tibble", "stringr", "data.table", "zoo", "igraph",
                  "openxlsx", "writexl", "RColorBrewer")

core_results <- sapply(core_packages, test_package)

# Bioconductor packages
cat("\nBioconductor Packages:\n")
cat("----------------------\n")
bioc_packages <- c("DESeq2", "edgeR", "sva", "biomaRt", "GenomicRanges",
                   "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm39.knownGene")

bioc_results <- sapply(bioc_packages, test_package)

# Special packages
cat("\nSpecial Packages:\n")
cat("-----------------\n")
special_results <- c(
  test_package("textshape"),
  test_package("plotgardener"),
  test_package("miqtl", "gkeele/miqtl")
)

# Summary
cat("\n====================================\n")
cat("Summary:\n")
cat("====================================\n")

total_core <- sum(core_results)
total_bioc <- sum(bioc_results)
total_special <- sum(special_results)
total_packages <- length(core_results) + length(bioc_results) + length(special_results)
total_loaded <- total_core + total_bioc + total_special

cat(sprintf("Core packages:        %d/%d loaded\n", total_core, length(core_results)))
cat(sprintf("Bioconductor packages: %d/%d loaded\n", total_bioc, length(bioc_results)))
cat(sprintf("Special packages:     %d/%d loaded\n", total_special, length(special_results)))
cat(sprintf("\nTotal:                %d/%d packages available\n", total_loaded, total_packages))

if (total_loaded == total_packages) {
  cat("\n✓ All packages loaded successfully!\n")
  cat("The environment is ready for CC-miQTL analysis.\n")
} else {
  cat("\n⚠ Some packages are missing.\n")
  cat("Run the setup script or install missing packages manually.\n")
}

# Test critical file access
cat("\n====================================\n")
cat("Testing File Access:\n")
cat("====================================\n")

genome_cache <- "/proj/raulab/projects/cc_gwas/data/raw/genomes/Corrected_Haplotypes_Cache_CC_83024"
if (file.exists(genome_cache)) {
  cat(sprintf("✓ Genome cache accessible: %s\n", genome_cache))
} else {
  cat(sprintf("✗ Genome cache NOT found: %s\n", genome_cache))
}

helper_script <- "/proj/raulab/projects/cc_gwas/scripts/extras/scan_h2lmm.R"
if (file.exists(helper_script)) {
  cat(sprintf("✓ Helper script accessible: %s\n", helper_script))
} else {
  cat(sprintf("✗ Helper script NOT found: %s\n", helper_script))
}

cat("\n")