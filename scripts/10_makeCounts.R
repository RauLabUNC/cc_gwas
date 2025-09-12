# =============================================================================
# Make a counts matrix from the raw STAR outputs
#
# Description:
# This script merges the barcode and feature information from STAR for 
# five plates, the n
# =============================================================================

# --- Load Libraries ---
library(tidyr)
library(dplyr)
library(stringr)
library(edgeR)
library(Matrix)
library(data.table)
library(tidyverse)
library(biomaRt)
library(optparse)

option_list <- list(
  make_option(c("--input_genes_in_trait_loci"), type = "character", help = "Path to input table of genes in trait loci"),
  make_option(c("--output_raw_counts"), type = "character", help = "Path to output csv of raw counts"),
  make_option(c("--output_counts_with_meta"), type = "character", help = "Path to output csv of counts with metadata")#,
  #make_option(c("--output_counts_for_eqtl"), type = "character", help = "Path to output RDS files of counts for eQTL")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)

print(opt)

genes_in_trait_loci <- readRDS(opt$options$input_genes_in_trait_loci)
str(genes_in_trait_loci)

# Function to read a plate's data
read_plate_data <- function(plate_number) {
  plate_dir <- paste0("data/processed/star_aligned/CC_Plate_", plate_number)
  # Read the sparse matrix
  mat <- readMM(paste0(plate_dir, "/Solo.out/Gene/raw/matrix.mtx"))
  
  # Read features (genes) and barcodes
  features <- fread(paste0(plate_dir, "/Solo.out/Gene/raw/features.tsv"), 
                    header=FALSE, stringsAsFactors=FALSE)
  barcodes <- fread(paste0(plate_dir, "/Solo.out/Gene/raw/barcodes.tsv"), 
                    header=FALSE, stringsAsFactors=FALSE)
  
  # Add plate identifier to barcodes
  barcodes$V1 <- paste0(plate_number, "_", barcodes$V1)
  
  # Convert to dense matrix and add row/column names
  mat <- as.matrix(mat)
  rownames(mat) <- features$V1
  colnames(mat) <- barcodes$V1
  
  return(mat)
}

# Read each plate's data. This is the count matrix, colnames = the plate # _ barcode, rownmes = transcripts
counts <- lapply(1:5, read_plate_data) %>%
  do.call(cbind, .) |> 
  as.data.frame()

# Load the Ensembl dataset for mouse genes
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract the Ensembl gene IDs from the row names
ensembl_ids <- rownames(counts)

# Query biomaRt to get the gene symbols for these Ensembl IDs
genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
               filters = 'ensembl_gene_id', 
               values = ensembl_ids, 
               mart = mart)

# Join the gene symbols back to the counts data frame
# First, make sure the Ensembl IDs are a column in the counts data frame for merging
counts$ensembl_gene_id <- rownames(counts)

# Merge to add gene symbols
counts_with_symbols <- left_join(counts, genes, by = 'ensembl_gene_id')


# First, identify numeric columns once
numeric_cols <- names(counts_with_symbols)[sapply(counts_with_symbols, is.numeric)]

# Convert to data.table
dt <- as.data.table(counts_with_symbols)

# Perform the summarization
summarized_counts <- dt[!is.na(external_gene_name),
                        lapply(.SD, sum, na.rm = TRUE),
                        by = external_gene_name,
                        .SDcols = numeric_cols]

# Replace Ensembl IDs with gene symbols where available, keeping Ensembl ID where not
rownames(summarized_counts) <- summarized_counts$external_gene_name

summarized_counts <- summarized_counts |> dplyr::select(-external_gene_name)

output_dir <- dirname(opt$options$output_raw_counts)
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
write.csv(summarized_counts, opt$options$output_raw_counts, row.names = T)

###Anh's part, clean up count matrix more### 
# srun -t 5:00:00 -p interact -n 1 --cpus-per-task=1 --mem=32g --pty /bin/bash
# summarized_counts <- read.csv("data/processed/joinLoci/relational_tables/raw_counts.csv", row.names = 1, check.names = FALSE)
# Import Brian's cleaned up sample info with ID, plate number, barcode, etc
sampleInfoBasic <- read.csv("data/processed/phenotypes/sampleInfoBasic.csv")
print("sampleInfoBasic str:")
str(sampleInfoBasic)
#Make a column that matches count colnames
sampleInfo = sampleInfoBasic %>%
  mutate(CountRep = paste(PlateNumber, AlitheaBarcode, sep="_"))

# Total reads per replicate (column) without transposing the whole matrix
total_reads <- colSums(summarized_counts)              # names = CountRep

rep_depth_df <- tibble(
  CountRep   = names(total_reads),
  TotalReads = as.numeric(total_reads)
)

print("rep_depth_df str:")
str(rep_depth_df)

# Join to sample metadata and choose the top-depth replicate per SampleID
rep_choice <- sampleInfo %>%
  inner_join(rep_depth_df, by = "CountRep") %>%
  group_by(SampleID) %>%
  slice_max(order_by = TotalReads, n = 1, with_ties = FALSE) %>%
  ungroup()

print("rep_choice str:")
str(rep_choice)

print("summarized_counts str:")
str(summarized_counts)

# Subset count matrix to the chosen replicates
keep <- colnames(summarized_counts) %in% rep_choice$CountRep # sanity check, should be all TRUE
filtered_counts <- as.data.frame(summarized_counts)[, keep]

# Reformat to sample rows + metadata
print("filtered_counts str:")
str(filtered_counts)
counts_long <- as.data.frame(t(filtered_counts))
counts_long$CountRep <- rownames(counts_long)

counts_info <- rep_choice %>%
  dplyr::select(SampleID, Strain, Sex, Drug, CountRep, Batch, TotalReads) %>%
  inner_join(counts_long, by = "CountRep") %>%
  relocate(SampleID, Strain, Sex, Drug, CountRep, Batch, TotalReads)

print("counts_info str:")
str(counts_info)
# Write output
write.csv(
  counts_info,
  opt$options$output_counts_with_meta,
  row.names = T
)
