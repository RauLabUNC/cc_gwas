# =============================================================================
# Consolidate Unique Position IDs from QTL Analyses
#
# Description:
# This script consolidates all unique positional IDs from trait and expression
# QTL tables, creating a unified position reference table with parsed
# chromosome coordinates. This is to prevent wasting time querying
# the same positions multiple times.
# =============================================================================

# --- Load Libraries ---
library(data.table)

# --- Read Position IDs from QTL Tables ---
trait_pos <- fread(
  "data/processed/joinLoci/relational_tables/traitLoci.csv",
  select = "Position_ID"  # Read only the needed column
)

eqtl_pos <- fread(
  "data/processed/joinLoci/relational_tables/expLoci.csv",
  select = "Position_ID"
)

# --- Combine and Parse Position IDs ---
# Union and deduplicate positions
pos <- unique(rbindlist(list(trait_pos, eqtl_pos)))

# Parse Position_ID format: "chr:start-end"
pos[, c("chr", "start_bp", "end_bp") := 
      tstrsplit(Position_ID, "[:-]", type.convert = TRUE, keep = 1:3)]

# Calculate locus length
pos[, length_bp := end_bp - start_bp + 1L]

# Sort by genomic position
setorder(pos, chr, start_bp)

# --- Format and Save Output ---
# Rename column for consistency
setnames(pos, "Position_ID", "pos_id")

# Ensure output directory exists
output_dir <- "data/processed/joinLoci/relational_tables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save position reference table
fwrite(pos, file.path(output_dir, "pos.csv"))

# Print summary statistics
cat("Unique positions found:", nrow(pos), "\n")
cat("Chromosomes covered:", length(unique(pos$chr)), "\n")
