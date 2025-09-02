# =============================================================================
# Split Gene Information into Species-Specific Tables
#
# Description:
# This script splits the comprehensive gene_info table from getPathogenicity.R
# into separate species-specific tables for mouse genes and orthology 
# relationships.
# =============================================================================

# --- Load Libraries ---
library(data.table)

# --- Load Gene Information ---
gi <- fread("data/processed/joinLoci/relational_tables/gene_info.csv",
            select = c("mouse_ensembl_id", "mouse_gene_symbol",
                       "description", "gene_biotype",
                       "chromosome_name", "start_position", "end_position",
                       "human_ensembl_id", "human_gene_symbol"))

cat("Loaded gene info for", nrow(gi), "genes\n")

# --- Create Mouse Gene Table ---
genes_mouse <- unique(
  gi[, .(mouse_ensembl_id,
         mouse_gene_symbol,
         description,
         gene_biotype,
         chr = chromosome_name,
         start_bp = start_position,
         end_bp = end_position)]
)

# Sort by genomic position
setorder(genes_mouse, chr, start_bp)
cat("Mouse genes table:", nrow(genes_mouse), "unique genes\n")

# --- Create Orthology Table ---
orthology <- unique(
  gi[, .(mouse_ensembl_id,
         human_ensembl_id,
         human_gene_symbol)]
)

cat("Orthology table:", nrow(orthology), "mappings\n")

# --- Save Output Tables ---
# Ensure output directory exists
output_dir <- "data/processed/joinLoci/relational_tables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save tables
fwrite(genes_mouse, file.path(output_dir, "genes_mouse.csv"))
fwrite(orthology, file.path(output_dir, "orthology.csv"))

cat("\nFiles saved successfully:\n")
cat("  - genes_mouse.csv\n")
cat("  - orthology.csv\n")