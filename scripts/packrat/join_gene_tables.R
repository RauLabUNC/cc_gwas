#' PackRat: Gene Table Joining Functions
#'
#' Functions for joining relational gene annotation tables (orthology, variants,
#' phenotypes, expression, eQTL) into unified gene summaries for locus packets.
#'
#' Based on the relational data structure shown in Figure B of Gural et al.
#'
#' @author Brian Gural, Anh Luu, Todd Kimball, Christoph Rau
#' @references Gural et al. (in preparation)

library(data.table)
library(dplyr)

# =============================================================================
# Core Gene Table Assembly
# =============================================================================

#' Build master gene table for a locus
#'
#' Joins genes in a locus with optional annotation tables (orthology, expression,
#' variants, eQTL, phenotypes). Designed to be flexible - accepts any combination
#' of optional data sources.
#'
#' @param genes_in_locus Data frame/table with gene info (must have gene ID column)
#' @param gene_id_col Character, name of gene ID column (e.g., "mouse_ensembl_id")
#' @param gene_symbol_col Character, name of gene symbol column (e.g., "mouse_gene_symbol")
#' @param orthology Optional data frame with orthology relationships
#' @param variants Optional data frame with coding variants per gene
#' @param expression Optional data frame with expression values (e.g., CPM, TPM)
#' @param eqtl Optional data frame with eQTL status per gene
#' @param phenotypes Optional data frame with phenotype annotations
#' @param custom_data Optional list of additional data frames to join
#'
#' @return Data table with merged gene information
#'
#' @details
#' This function performs sequential left joins of all provided data sources
#' onto the base gene table. Each optional table should contain a column
#' matching `gene_id_col` or `gene_symbol_col` for joining.
#'
#' The function is intentionally flexible to accommodate different:
#' - Model organisms (mouse, rat, fly, human)
#' - Gene identifiers (Ensembl IDs, gene symbols, Entrez IDs)
#' - Data sources (user-specific expression, eQTL, variant calling)
#'
#' @examples
#' # Minimal example - genes only
#' genes <- data.frame(
#'   mouse_ensembl_id = c("ENSMUSG001", "ENSMUSG002"),
#'   mouse_gene_symbol = c("Abcb10", "Acsl5"),
#'   chr = c(8, 8),
#'   start_bp = c(1000000, 2000000),
#'   end_bp = c(1050000, 2050000)
#' )
#' master_table <- build_gene_table(genes, gene_id_col = "mouse_ensembl_id")
#'
#' # With expression data
#' expr <- data.frame(
#'   mouse_gene_symbol = c("Abcb10", "Acsl5"),
#'   avg_cpm_ctrl = c(150.2, 45.8)
#' )
#' master_table <- build_gene_table(genes,
#'                                  gene_id_col = "mouse_ensembl_id",
#'                                  gene_symbol_col = "mouse_gene_symbol",
#'                                  expression = expr)
#'
#' @export
build_gene_table <- function(genes_in_locus,
                              gene_id_col = "mouse_ensembl_id",
                              gene_symbol_col = "mouse_gene_symbol",
                              orthology = NULL,
                              variants = NULL,
                              expression = NULL,
                              eqtl = NULL,
                              phenotypes = NULL,
                              custom_data = NULL) {

  # Convert to data.table for efficient joining
  setDT(genes_in_locus)
  master_table <- copy(genes_in_locus)

  # Join orthology (if provided)
  if (!is.null(orthology)) {
    setDT(orthology)
    # Try joining by gene ID first, then gene symbol
    if (gene_id_col %in% names(orthology)) {
      master_table <- merge(master_table, orthology,
                           by = gene_id_col, all.x = TRUE)
    } else if (gene_symbol_col %in% names(orthology)) {
      master_table <- merge(master_table, orthology,
                           by = gene_symbol_col, all.x = TRUE)
    } else {
      warning("Orthology table does not contain matching gene ID column. Skipping.")
    }
  }

  # Join variants (if provided)
  if (!is.null(variants)) {
    setDT(variants)
    if (gene_symbol_col %in% names(variants)) {
      master_table <- merge(master_table, variants,
                           by = gene_symbol_col, all.x = TRUE)
    } else if (gene_id_col %in% names(variants)) {
      master_table <- merge(master_table, variants,
                           by = gene_id_col, all.x = TRUE)
    } else {
      warning("Variants table does not contain matching gene column. Skipping.")
    }
  }

  # Join expression (if provided)
  if (!is.null(expression)) {
    setDT(expression)
    if (gene_symbol_col %in% names(expression)) {
      master_table <- merge(master_table, expression,
                           by = gene_symbol_col, all.x = TRUE)
    } else if (gene_id_col %in% names(expression)) {
      master_table <- merge(master_table, expression,
                           by = gene_id_col, all.x = TRUE)
    } else {
      warning("Expression table does not contain matching gene column. Skipping.")
    }
  }

  # Join eQTL (if provided)
  if (!is.null(eqtl)) {
    setDT(eqtl)
    if (gene_symbol_col %in% names(eqtl)) {
      master_table <- merge(master_table, eqtl,
                           by = gene_symbol_col, all.x = TRUE)
    } else if (gene_id_col %in% names(eqtl)) {
      master_table <- merge(master_table, eqtl,
                           by = gene_id_col, all.x = TRUE)
    } else {
      warning("eQTL table does not contain matching gene column. Skipping.")
    }
  }

  # Join phenotypes (if provided)
  if (!is.null(phenotypes)) {
    setDT(phenotypes)
    if (gene_symbol_col %in% names(phenotypes)) {
      master_table <- merge(master_table, phenotypes,
                           by = gene_symbol_col, all.x = TRUE)
    } else if (gene_id_col %in% names(phenotypes)) {
      master_table <- merge(master_table, phenotypes,
                           by = gene_id_col, all.x = TRUE)
    } else {
      warning("Phenotypes table does not contain matching gene column. Skipping.")
    }
  }

  # Join custom data (if provided)
  if (!is.null(custom_data)) {
    if (!is.list(custom_data)) {
      warning("custom_data must be a list of data frames. Skipping.")
    } else {
      for (i in seq_along(custom_data)) {
        custom_df <- custom_data[[i]]
        setDT(custom_df)

        # Try joining by gene symbol, then gene ID
        if (gene_symbol_col %in% names(custom_df)) {
          master_table <- merge(master_table, custom_df,
                               by = gene_symbol_col, all.x = TRUE)
        } else if (gene_id_col %in% names(custom_df)) {
          master_table <- merge(master_table, custom_df,
                               by = gene_id_col, all.x = TRUE)
        } else {
          warning(paste0("Custom data frame ", i,
                        " does not contain matching gene column. Skipping."))
        }
      }
    }
  }

  return(master_table)
}


# =============================================================================
# Helper Functions for Specific Joins
# =============================================================================

#' Get genes within a genomic region
#'
#' @param chr Character or numeric, chromosome identifier
#' @param start Numeric, region start position (bp)
#' @param end Numeric, region end position (bp)
#' @param gene_data Data table with columns: chr, start_bp, end_bp
#' @param chr_col Name of chromosome column (default "chr")
#' @param start_col Name of start position column (default "start_bp")
#' @param end_col Name of end position column (default "end_bp")
#'
#' @return Data table of genes overlapping the region
#'
#' @examples
#' genes <- data.frame(
#'   gene = c("Abcb10", "Acsl5", "Mir1967"),
#'   chr = c(8, 8, 8),
#'   start_bp = c(1000000, 2000000, 3000000),
#'   end_bp = c(1050000, 2050000, 3001000)
#' )
#' get_genes_in_region(chr = 8, start = 900000, end = 2100000, gene_data = genes)
#' # Returns Abcb10 and Acsl5
#'
#' @export
get_genes_in_region <- function(chr, start, end, gene_data,
                                 chr_col = "chr",
                                 start_col = "start_bp",
                                 end_col = "end_bp") {
  setDT(gene_data)
  gene_data[get(chr_col) == chr &
            get(start_col) < end &
            get(end_col) > start, ]
}


#' Add locus membership columns for multiple loci
#'
#' For a gene table, add binary (Yes/No) columns indicating whether each gene
#' is present in each locus. Useful for multi-locus/cluster packets.
#'
#' @param gene_table Data table with gene coordinates
#' @param loci_df Data frame with locus definitions (chr, start, end, id)
#' @param chr_col Name of chromosome column in gene_table
#' @param start_col Name of start position column in gene_table
#' @param end_col Name of end position column in gene_table
#' @param locus_id_col Name of column to use for locus identifiers in loci_df
#'
#' @return gene_table with added "In_[locus_id]" columns (Yes/No)
#'
#' @details
#' This creates columns like "In_LV_chr8_28.3Mb" to indicate which genes
#' fall within each locus in a cluster of overlapping QTL.
#'
#' @export
add_locus_membership <- function(gene_table, loci_df,
                                  chr_col = "chr",
                                  start_col = "start_bp",
                                  end_col = "end_bp",
                                  locus_id_col = "locus_id") {

  setDT(gene_table)
  setDT(loci_df)

  # Ensure loci have ID column
  if (!locus_id_col %in% names(loci_df)) {
    stop(paste0("loci_df must contain a '", locus_id_col, "' column"))
  }

  # For each locus, add a membership column
  for (i in 1:nrow(loci_df)) {
    locus <- loci_df[i, ]
    locus_col_name <- paste0("In_", locus[[locus_id_col]])

    # Determine locus boundaries
    locus_chr <- locus$chr
    locus_start <- locus$start
    locus_end <- locus$end

    # Initialize column to "No"
    gene_table[, (locus_col_name) := "No"]

    # Set to "Yes" for genes in this locus
    gene_table[get(chr_col) == locus_chr &
               get(start_col) < locus_end &
               get(end_col) > locus_start,
               (locus_col_name) := "Yes"]
  }

  return(gene_table)
}


#' Summarize variants per gene
#'
#' Collapses variant-level data into gene-level summaries (e.g., concatenate
#' all variant types per gene).
#'
#' @param variants_df Data frame with columns: gene (symbol or ID), variant info
#' @param gene_col Name of gene column in variants_df
#' @param variant_col Name of column with variant descriptions
#' @param sep Separator for concatenating multiple variants (default "; ")
#'
#' @return Data table with one row per gene and concatenated variant info
#'
#' @examples
#' variants <- data.frame(
#'   gene = c("Abcb10", "Abcb10", "Acsl5"),
#'   mutation = c("Missense_WSB", "Missense_CAST", "StopLost_PWK")
#' )
#' summarize_variants(variants, gene_col = "gene", variant_col = "mutation")
#'
#' @export
summarize_variants <- function(variants_df,
                                gene_col = "gene",
                                variant_col = "mutation",
                                sep = "; ") {
  setDT(variants_df)
  variants_summary <- variants_df[, .(
    variants = paste(get(variant_col), collapse = sep)
  ), by = gene_col]

  # Rename back to original column name
  setnames(variants_summary, "gene_col", gene_col)
  return(variants_summary)
}


# =============================================================================
# Column Ordering and Formatting
# =============================================================================

#' Reorder columns in gene table for readability
#'
#' Puts essential columns first, locus membership columns next, then remaining.
#'
#' @param gene_table Data table to reorder
#' @param essential_cols Character vector of column names to put first
#' @param locus_col_pattern Regex pattern to identify locus membership columns (default "^In_")
#'
#' @return Data table with reordered columns
#'
#' @export
reorder_gene_columns <- function(gene_table,
                                  essential_cols = c("mouse_gene_symbol",
                                                     "mouse_ensembl_id",
                                                     "chr", "start_bp", "end_bp"),
                                  locus_col_pattern = "^In_") {

  setDT(gene_table)

  # Identify locus membership columns
  locus_cols <- grep(locus_col_pattern, names(gene_table), value = TRUE)

  # Identify remaining columns
  remaining_cols <- setdiff(names(gene_table), c(essential_cols, locus_cols))

  # Only include essential_cols that actually exist in the table
  essential_cols <- intersect(essential_cols, names(gene_table))

  # Create final column order
  col_order <- c(essential_cols, locus_cols, remaining_cols)

  return(gene_table[, ..col_order])
}
