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
library(optparse)

option_list <- list(
  make_option(c("--input_gene_info"), type = "character", help = "Path to input table of genes in trait loci"),
  make_option(c("--output_mouse_genes"), type = "character", help = "Path to output csv of mouse genes"),
  make_option(c("--output_orthology"), type = "character", help = "Path to output csv of orthology relationships")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)

print(opt)
# Old input path
# "data/processed/joinLoci/relational_tables/gene_info.csv"
# new path

# Old outs
# "data/processed/joinLoci/relational_tables/genes_mouse.csv"
# "data/processed/joinLoci/relational_tables/orthology.csv"
# --- Load Gene Information ---
gi <- fread(opt$options$input_gene_info,
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
output_dir <- dirname(opt$options$output_mouse_genes)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save tables
fwrite(genes_mouse, opt$options$output_mouse_genes)
fwrite(orthology, opt$options$output_orthology)

cat("\nFiles saved successfully:\n")
cat("  - genes_mouse.csv\n")
cat("  - orthology.csv\n")