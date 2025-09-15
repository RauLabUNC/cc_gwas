# --- Load Libraries ---
library(miqtl)
library(tidyverse)
library(optparse)
library(evir)

# --- Define Command-Line Arguments ---
option_list <- list(
  make_option(c("--output_threshold"), type = "character", help = "Path to final output threshold RDS file"),
  make_option(c("--output_max_stats"), type = "character", help = "Path to final output combined max statistics RDS file")
)

# --- Parse Arguments ---
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)
chunk_files <- opt$args

str(opt)

# --- Load and Combine Chunks ---
cat(sprintf("Loading and combining max statistics from %d chunks...\n", length(chunk_files)))
all_max_stats_chunks <- lapply(chunk_files, readRDS)


combined_stats <- list(
  LOD = do.call(c, lapply(all_max_stats_chunks, `[[`, "LOD")),
  p.value = do.call(c, lapply(all_max_stats_chunks, `[[`, "p.value"))
)

# Create a mock object for the threshold function.
mock_threshold_scans <- list(
  max.statistics = combined_stats
)

# --- Calculate Final Threshold ---
cat("Calculating GEV thresholds from combined statistics...\n")
permute_threshold <- get.gev.thresholds(
  threshold.scans = mock_threshold_scans,
  percentile = 0.85,
  use.lod = TRUE
)

print(opt$options$output_max_stats)

# --- Save Final Outputs ---
output_dir <- dirname(opt$options$output_max_stats)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

saveRDS(permute_threshold, opt$options$output_threshold)
saveRDS(combined_stats, opt$options$output_max_stats)