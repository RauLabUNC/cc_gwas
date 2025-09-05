# --- Load Libraries ---
library(miqtl)
library(tidyverse)
library(optparse)

# --- Define Command-Line Arguments ---
option_list <- list(
  make_option(c("--input_perms"), type = "character", help = "Path to input permuted phenotypes RDS file"),
  make_option(c("--input_pheno"), type = "character", help = "Path to input processed phenotype file"),
  make_option(c("--output_scan_chunk"), type = "character", help = "Path to output permuted scan chunk RDS file"),
  make_option(c("--drug"), type = "character", help = "Drug treatment (e.g., Ctrl, Iso)"),
  make_option(c("--qtl_trait"), type = "character", help = "QTL trait to analyze"),
  make_option(c("--chunk_index"), type = "integer", help = "Index of the permutation chunk to process"),
  make_option(c("--chunk_size"), type = "integer", help = "Number of permutations per chunk"),
  make_option(c("--mode"), type = "character", default = "full", help = "Run mode: 'full' or 'test'")
)

# Parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# --- Load Input Data ---
permuted_phenotype <- readRDS(opt$input_perms)
phenotypes <- read.csv(opt$input_pheno)
phenotypes <- phenotypes |> filter(!is.na(opt$qtl_trait) & Drug == opt$drug)

# --- Subset Permutations ---
total_perms <- ncol(permuted_phenotype$y.matrix)
start_idx <- (opt$chunk_index - 1) * opt$chunk_size + 1
end_idx <- min(opt$chunk_index * opt$chunk_size, total_perms)

if (start_idx > total_perms) {
  stop("Chunk index is out of bounds.")
}

cat(sprintf("Processing permutation chunk %d: samples %d to %d\n", opt$chunk_index, start_idx, end_idx))

# Subset the sim.threshold.object
permuted_phenotype_chunk <- permuted_phenotype
permuted_phenotype_chunk$y.matrix <- permuted_phenotype_chunk$y.matrix[, start_idx:end_idx, drop = FALSE]

message("Permuted phenotype chunk structure:")
message(permuted_phenotype_chunk)
# --- Genome Scan Setup ---
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"
chr_to_scan <- if (opt$mode == "test") 1 else "all"

# Run scans on the chunk of permuted phenotypes
permuted_scans_chunk <- run.threshold.scans(
  sim.threshold.object = permuted_phenotype_chunk,
  keep.full.scans = FALSE,
  genomecache = genomecache,
  data = phenotypes,
  use.multi.impute = FALSE,
  scan.seed = opt$chunk_index, # Use chunk index for a unique seed
  chr = chr_to_scan
)

print("Permuted scans chunk structure:")
message(str(permuted_scans_chunk))
# Extract just the max statistics to save space
max_stats_chunk <- permuted_scans_chunk$max.statistics
message("Max stats chunk structure:")
message(str(max_stats_chunk))

# --- Save Output ---
output_dir <- dirname(opt$output_scan_chunk)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(max_stats_chunk, opt$output_scan_chunk)