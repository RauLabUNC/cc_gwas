#' PackRat: Summary Report Generation
#'
#' Functions for creating README files and text summaries for locus packets.
#' These provide high-level overviews for rapid assessment of candidate genes.
#'
#' @author Brian Gural, Anh Luu, Todd Kimball, Christoph Rau
#' @references Gural et al. (in preparation)

library(dplyr)

# =============================================================================
# Main README Generation
# =============================================================================

#' Generate README summary for a locus packet
#'
#' Creates a formatted text file summarizing key information about a locus
#' and its candidate genes. Customizable for different study types.
#'
#' @param locus_info List or data frame row with locus information
#'   Should contain: chr, start, end, trait (optional), peak (optional)
#' @param gene_summary Data frame with gene summary statistics
#' @param output_file Character, path to output README file
#' @param phenotype_keywords Character vector of keywords used for filtering (optional)
#' @param custom_sections List of custom text sections to include (optional)
#'   Each element should be named (section title) and contain text content
#' @param genome_build Character, genome build identifier (e.g., "mm39", "GRCh38")
#'
#' @return Invisibly returns the README text
#'
#' @details
#' The README provides:
#' - Locus coordinates and associated traits
#' - Gene counts (total, filtered, with supporting evidence)
#' - Lists of top candidate genes
#' - Optional custom sections
#'
#' @examples
#' locus <- list(chr = 8, start = 28300000, end = 41100000,
#'              trait = "LV_Mass", peak = 34500000)
#' genes <- data.frame(gene = c("Abcb10", "Acsl5"), has_eqtl = c("Yes", "No"))
#' create_locus_readme(locus, genes, "README_locus.txt")
#'
#' @export
create_locus_readme <- function(locus_info,
                                gene_summary = NULL,
                                output_file = "README_summary.txt",
                                phenotype_keywords = NULL,
                                custom_sections = NULL,
                                genome_build = "mm39") {

  # Extract locus information
  chr <- locus_info$chr
  start_mb <- locus_info$start / 1e6
  end_mb <- locus_info$end / 1e6
  trait <- if (!is.null(locus_info$trait)) locus_info$trait else "Unknown"
  peak_mb <- if (!is.null(locus_info$peak)) locus_info$peak / 1e6 else NA

  # Build summary text
  summary_lines <- c(
    "# LOCUS SUMMARY",
    paste0("Genomic Region: chr", chr, ":", start_mb, "-", end_mb, " Mb (", genome_build, ")"),
    if (!is.na(peak_mb)) paste0("Peak Position: ", round(peak_mb, 2), " Mb") else NULL,
    if (!is.null(trait)) paste0("Associated Trait: ", trait) else NULL,
    ""
  )

  # Gene statistics (if provided)
  if (!is.null(gene_summary)) {
    total_genes <- nrow(gene_summary)
    summary_lines <- c(summary_lines,
                      "# GENE STATISTICS",
                      paste0("Total genes in region: ", total_genes))

    # Count genes with different types of evidence (if columns exist)
    if ("has_orthologue" %in% names(gene_summary)) {
      n_ortho <- sum(gene_summary$has_orthologue == "Yes", na.rm = TRUE)
      summary_lines <- c(summary_lines,
                        paste0("  - With human orthologue: ", n_ortho,
                              " (", round(100 * n_ortho / total_genes, 1), "%)"))
    }

    if ("protein_coding" %in% names(gene_summary)) {
      n_protein <- sum(gene_summary$protein_coding == "Yes", na.rm = TRUE)
      summary_lines <- c(summary_lines,
                        paste0("  - Protein coding: ", n_protein,
                              " (", round(100 * n_protein / total_genes, 1), "%)"))
    }

    if ("has_eqtl" %in% names(gene_summary)) {
      n_eqtl <- sum(gene_summary$has_eqtl == "Yes", na.rm = TRUE)
      summary_lines <- c(summary_lines,
                        paste0("  - With cis-eQTL: ", n_eqtl,
                              " (", round(100 * n_eqtl / total_genes, 1), "%)"))
    }

    if ("has_variant" %in% names(gene_summary)) {
      n_variant <- sum(gene_summary$has_variant == "Yes", na.rm = TRUE)
      summary_lines <- c(summary_lines,
                        paste0("  - With coding variant: ", n_variant,
                              " (", round(100 * n_variant / total_genes, 1), "%)"))
    }

    if ("differentially_expressed" %in% names(gene_summary)) {
      n_de <- sum(gene_summary$differentially_expressed == "Yes", na.rm = TRUE)
      summary_lines <- c(summary_lines,
                        paste0("  - Differentially expressed: ", n_de,
                              " (", round(100 * n_de / total_genes, 1), "%)"))
    }

    summary_lines <- c(summary_lines, "")
  }

  # Add phenotype keyword info (if provided)
  if (!is.null(phenotype_keywords)) {
    summary_lines <- c(summary_lines,
                      "# FILTERING CRITERIA",
                      paste0("Phenotype keywords used: ",
                            paste(head(phenotype_keywords, 10), collapse = ", "),
                            if (length(phenotype_keywords) > 10) "..." else ""),
                      "")
  }

  # Add custom sections
  if (!is.null(custom_sections)) {
    for (section_name in names(custom_sections)) {
      summary_lines <- c(summary_lines,
                        paste0("# ", toupper(section_name)),
                        custom_sections[[section_name]],
                        "")
    }
  }

  # Add footer
  summary_lines <- c(summary_lines,
                    "# NOTES",
                    paste0("Generated: ", Sys.time()),
                    "See accompanying Excel file for detailed gene annotations.")

  # Combine into single text
  summary_text <- paste(summary_lines, collapse = "\n")

  # Write to file
  message("Creating README: ", output_file)
  writeLines(summary_text, output_file)

  invisible(summary_text)
}


# =============================================================================
# Specialized Summary Functions
# =============================================================================

#' Create phenotype-focused README (e.g., cardiac genes)
#'
#' Generates a README highlighting genes with specific phenotype associations.
#'
#' @param locus_info List with locus information
#' @param gene_table Data frame with gene annotations
#' @param phenotype_keywords Character vector of phenotype keywords
#' @param output_file Character, output file path
#' @param gene_symbol_col Name of gene symbol column
#' @param human_disease_col Name of human disease summary column
#' @param mouse_pheno_col Name of mouse phenotype summary column
#' @param ... Additional arguments passed to create_locus_readme()
#'
#' @return Invisibly returns README text
#'
#' @details
#' This creates a README with a section specifically highlighting genes that
#' match phenotype keywords in either human disease or mouse phenotype databases.
#' Useful for study-specific prioritization.
#'
#' @export
create_phenotype_readme <- function(locus_info,
                                    gene_table,
                                    phenotype_keywords,
                                    output_file = "README_summary.txt",
                                    gene_symbol_col = "mouse_gene_symbol",
                                    human_disease_col = "Human Disease Summary",
                                    mouse_pheno_col = "Mouse Phenotype Summary (MGI)",
                                    ...) {

  # Source the phenotype extraction functions
  source("scripts/packrat/phenotype_extract.R")

  # Build keyword pattern
  pattern <- paste0("\\b(", paste(phenotype_keywords, collapse = "|"), ")\\b")

  # Filter genes with matching phenotypes
  human_matches <- gene_table %>%
    filter(grepl(pattern, get(human_disease_col), ignore.case = TRUE))

  mouse_matches <- gene_table %>%
    filter(grepl(pattern, get(mouse_pheno_col), ignore.case = TRUE))

  # Create summary sections
  custom_sections <- list()

  if (nrow(human_matches) > 0 || nrow(mouse_matches) > 0) {
    pheno_summary <- c(
      paste0("Genes with matching phenotypes: ",
            length(unique(c(human_matches[[gene_symbol_col]],
                          mouse_matches[[gene_symbol_col]])))),
      paste0("  - Human disease associations: ", nrow(human_matches)),
      paste0("  - Mouse phenotypes: ", nrow(mouse_matches))
    )

    # List top genes
    if (nrow(human_matches) > 0) {
      top_human <- head(human_matches[[gene_symbol_col]], 10)
      pheno_summary <- c(pheno_summary,
                        "",
                        "Top genes with human disease associations:",
                        paste0("  ", top_human))
    }

    if (nrow(mouse_matches) > 0) {
      top_mouse <- head(mouse_matches[[gene_symbol_col]], 10)
      pheno_summary <- c(pheno_summary,
                        "",
                        "Top genes with mouse phenotypes:",
                        paste0("  ", top_mouse))
    }

    custom_sections[["PHENOTYPE MATCHES"]] <- paste(pheno_summary, collapse = "\n")
  }

  # Generate README with custom sections
  create_locus_readme(
    locus_info = locus_info,
    gene_summary = gene_table,
    output_file = output_file,
    phenotype_keywords = phenotype_keywords,
    custom_sections = custom_sections,
    ...
  )
}


# =============================================================================
# CSV Summary Tables
# =============================================================================

#' Create CSV summary of top candidate genes
#'
#' Exports a filtered table of top candidate genes based on supporting evidence.
#'
#' @param gene_table Data frame with gene annotations
#' @param output_file Character, path to output CSV
#' @param min_evidence Numeric, minimum number of supporting evidence types (default 2)
#' @param evidence_cols Character vector of column names representing evidence types
#'   (e.g., c("has_eqtl", "has_variant", "differentially_expressed"))
#'
#' @return Invisibly returns the filtered gene table
#'
#' @export
create_top_candidates_csv <- function(gene_table,
                                      output_file,
                                      min_evidence = 2,
                                      evidence_cols = c("has_eqtl",
                                                       "has_variant",
                                                       "differentially_expressed")) {

  # Count evidence for each gene
  evidence_present <- evidence_cols[evidence_cols %in% names(gene_table)]

  if (length(evidence_present) == 0) {
    warning("No evidence columns found in gene table. Returning all genes.")
    top_genes <- gene_table
  } else {
    gene_table$evidence_count <- rowSums(
      gene_table[, evidence_present] == "Yes", na.rm = TRUE
    )

    # Filter for genes with minimum evidence
    top_genes <- gene_table %>%
      filter(evidence_count >= min_evidence) %>%
      arrange(desc(evidence_count))
  }

  # Write to CSV
  message("Creating top candidates CSV: ", output_file)
  write.csv(top_genes, output_file, row.names = FALSE)

  invisible(top_genes)
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Format large numbers with commas
#'
#' @param x Numeric value
#' @return Character string with comma separators
#' @keywords internal
format_number <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}


#' Convert base pairs to megabases
#'
#' @param bp Numeric, position in base pairs
#' @param digits Integer, number of decimal places (default 2)
#' @return Numeric value in megabases
#' @keywords internal
bp_to_mb <- function(bp, digits = 2) {
  round(bp / 1e6, digits)
}


#' Create file-safe locus identifier
#'
#' @param chr Chromosome
#' @param start Start position (bp)
#' @param end End position (bp)
#' @param trait Optional trait name
#' @return Character string suitable for filenames
#'
#' @examples
#' create_locus_id(8, 28300000, 41100000, "LV_Mass")
#' # Returns "chr8_28-41Mb_LV_Mass"
#'
#' @export
create_locus_id <- function(chr, start, end, trait = NULL) {
  start_mb <- round(start / 1e6, 0)
  end_mb <- round(end / 1e6, 0)

  locus_id <- paste0("chr", chr, "_", start_mb, "-", end_mb, "Mb")

  if (!is.null(trait)) {
    # Clean trait name for filename
    clean_trait <- gsub("[^A-Za-z0-9_-]", "_", trait)
    locus_id <- paste0(locus_id, "_", clean_trait)
  }

  return(locus_id)
}
