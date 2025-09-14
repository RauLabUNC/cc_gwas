# =============================================================================
# Identify Genes Within QTL Loci
#
# Description:
# This script identifies genes that overlap with QTL loci using biomaRt
# gene annotations. It performs genomic overlap analysis and filters out
# predicted genes (Gm*, *Rik).
# =============================================================================

# --- Load Libraries ---
library(data.table)
library(dplyr)
library(tidyr)
library(biomaRt)
library(optparse)

option_list <- list(
  make_option(c("--input_pos_summary"), type = "character", help = "Paths to input RDS files of positions"),
  make_option(c("--output_genes"), type = "character", help = "Paths to output RDS files of genes")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)

print(opt)
# --- Read and Prepare Loci Table ---
pos_dt <- fread(
  opt$options$input_pos_summary,
  colClasses = list(
    character = "pos_id",
    character = "chr",
    integer = c("start_bp", "end_bp")
  )
) |> 
  drop_na()

# Rename columns for foverlaps compatibility
setnames(pos_dt,
         old = c("chr", "start_bp", "end_bp"),
         new = c("chromosome_name", "start", "end"))

# Set key for efficient overlap operations
setkey(pos_dt, chromosome_name, start, end)

# --- Fetch Gene Annotations from BioMart ---
# Connect to Ensembl
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

# Fetch gene coordinates and names
genes_df <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  mart = mart
)

# Filter out predicted genes and format
genes_dt <- as.data.table(genes_df)[
  !grepl("^Gm", external_gene_name) &  # Remove Gm* predicted genes
  !grepl("Rik$", external_gene_name) &  # Remove *Rik predicted genes
  external_gene_name != "",              # Remove empty gene names
  .(
    chromosome_name,
    start = start_position,
    end = end_position,
    ensembl_gene_id,
    gene_symbol = external_gene_name
  )
]

# Set key for overlap operations
setkey(genes_dt, chromosome_name, start, end)

cat("Genes fetched:", nrow(genes_dt), "\n")

# --- Perform Genomic Overlap Analysis ---
# Find genes that overlap with loci positions
overlap_dt <- foverlaps(
  genes_dt,
  pos_dt,
  nomatch = 0L  # Exclude non-overlapping entries
)

cat("Overlaps found:", nrow(overlap_dt), "\n")

# --- Create Gene Lists per Locus ---
# Aggregate genes for each locus
locus_list_dt <- overlap_dt[
  , .(
    Gene_Symbols = list(unique(gene_symbol)),
    Ensembl_IDs = list(unique(ensembl_gene_id))
  ), by = pos_id
]

# --- Save Output ---
# Ensure output directory exists
output_dir <- dirname(opt$options$output_genes)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save gene lists per locus
saveRDS(locus_list_dt, opt$options$output_genes)

# Print summary
cat("Loci with genes:", nrow(locus_list_dt), "\n")
cat("Average genes per locus:", 
    round(mean(sapply(locus_list_dt$Gene_Symbols, length)), 1), "\n")
