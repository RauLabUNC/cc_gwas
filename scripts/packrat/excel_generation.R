#' PackRat: Excel Workbook Generation
#'
#' Functions for creating multi-sheet Excel workbooks from gene tables.
#' Designed to be flexible and accept any combination of gene annotations.
#'
#' @author Brian Gural, Anh Luu, Todd Kimball, Christoph Rau
#' @references Gural et al. (in preparation)

library(openxlsx)
library(data.table)
library(dplyr)

# =============================================================================
# Main Excel Generation Function
# =============================================================================

#' Create multi-sheet gene information workbook
#'
#' Generates a standardized Excel file with gene annotations, optionally including
#' phenotype/disease associations, expression data, variants, etc.
#'
#' @param gene_table Data frame with master gene table (one row per gene)
#' @param output_file Character, path to output .xlsx file
#' @param main_sheet_name Character, name for the main gene summary sheet (default "AllGenes")
#' @param phenotype_tables Optional list of phenotype/disease data frames to add as separate sheets
#'   Each list element should be a data frame with a 'gene_id' or 'gene_symbol' column
#'   Names of list elements become sheet names (e.g., list(MousePhenotypes = df, HumanDiseases = df))
#' @param gene_id_col Character, name of gene ID column for joining phenotype tables
#' @param essential_cols Character vector of column names to place first (left-most)
#' @param highlight_cols Character vector of column names to apply conditional formatting
#'   (e.g., "cis_eqtl", "has_variant") - will highlight "Yes" values
#' @param freeze_panes Logical, freeze first row and optionally first column (default TRUE)
#' @param auto_filter Logical, add auto-filter to main sheet (default TRUE)
#' @param wrap_text Logical, wrap text in all cells (default TRUE)
#' @param auto_width Logical, auto-fit column widths (default TRUE)
#'
#' @return Invisibly returns the path to the created Excel file
#'
#' @details
#' This function creates a clean, filterable Excel workbook suitable for manual
#' gene curation. The main sheet contains the master gene table, and optional
#' sheets contain detailed phenotype/disease annotations.
#'
#' Conditional formatting highlights:
#' - "Yes" values in binary columns (light green)
#' - Can be customized with highlight_cols parameter
#'
#' @examples
#' # Basic usage - gene table only
#' genes <- data.frame(
#'   gene_symbol = c("Abcb10", "Acsl5"),
#'   chr = c(8, 8),
#'   expression_cpm = c(150, 45),
#'   has_eqtl = c("Yes", "No"),
#'   has_variant = c("Yes", "Yes")
#' )
#' create_gene_workbook(genes, "locus_chr8_genes.xlsx",
#'                      highlight_cols = c("has_eqtl", "has_variant"))
#'
#' # With phenotype sheets
#' mouse_pheno <- data.frame(
#'   gene_symbol = c("Abcb10", "Abcb10", "Acsl5"),
#'   phenotype = c("abnormal heart", "increased stress response", "normal"),
#'   pubmed_id = c("12345", "67890", "11111")
#' )
#' create_gene_workbook(genes, "locus_chr8_genes.xlsx",
#'                      phenotype_tables = list(MousePhenotypes = mouse_pheno),
#'                      gene_id_col = "gene_symbol")
#'
#' @export
create_gene_workbook <- function(gene_table,
                                  output_file,
                                  main_sheet_name = "AllGenes",
                                  phenotype_tables = NULL,
                                  gene_id_col = "mouse_ensembl_id",
                                  essential_cols = NULL,
                                  highlight_cols = NULL,
                                  freeze_panes = TRUE,
                                  auto_filter = TRUE,
                                  wrap_text = TRUE,
                                  auto_width = TRUE) {

  # Create workbook
  wb <- createWorkbook()

  # Ensure gene_table is data frame (not data.table, for openxlsx compatibility)
  gene_table <- as.data.frame(gene_table)

  # Reorder columns if essential_cols specified
  if (!is.null(essential_cols)) {
    essential_present <- intersect(essential_cols, names(gene_table))
    other_cols <- setdiff(names(gene_table), essential_present)
    gene_table <- gene_table[, c(essential_present, other_cols)]
  }

  # --- Add Main Gene Summary Sheet ---
  addWorksheet(wb, main_sheet_name)
  writeDataTable(wb, sheet = main_sheet_name, x = gene_table,
                 tableName = gsub("[^A-Za-z0-9]", "", main_sheet_name),
                 withFilter = auto_filter)

  # Apply conditional formatting to highlight columns
  if (!is.null(highlight_cols)) {
    for (col_name in highlight_cols) {
      if (col_name %in% names(gene_table)) {
        col_idx <- which(names(gene_table) == col_name)

        # Highlight "Yes" values in green
        conditionalFormatting(
          wb, sheet = main_sheet_name,
          cols = col_idx,
          rows = 2:(nrow(gene_table) + 1),  # +1 for header
          rule = '=="Yes"',
          style = createStyle(bgFill = "#90EE90")  # Light green
        )

        # Optionally highlight "No" values in light red
        conditionalFormatting(
          wb, sheet = main_sheet_name,
          cols = col_idx,
          rows = 2:(nrow(gene_table) + 1),
          rule = '=="No"',
          style = createStyle(bgFill = "#FFB6C1")  # Light red
        )
      }
    }
  }

  # Apply text wrapping
  if (wrap_text) {
    wrap_style <- createStyle(wrapText = TRUE)
    addStyle(wb, sheet = main_sheet_name, style = wrap_style,
             rows = 1:(nrow(gene_table) + 1),
             cols = 1:ncol(gene_table),
             gridExpand = TRUE)
  }

  # Auto-fit column widths
  if (auto_width) {
    setColWidths(wb, sheet = main_sheet_name,
                 cols = 1:ncol(gene_table),
                 widths = "auto")
  }

  # Freeze panes
  if (freeze_panes) {
    freezePane(wb, sheet = main_sheet_name, firstActiveRow = 2)
  }

  # --- Add Phenotype/Disease Sheets ---
  if (!is.null(phenotype_tables)) {
    if (!is.list(phenotype_tables)) {
      warning("phenotype_tables must be a named list. Skipping.")
    } else {
      for (sheet_name in names(phenotype_tables)) {
        pheno_df <- as.data.frame(phenotype_tables[[sheet_name]])

        if (nrow(pheno_df) > 0) {
          # Truncate sheet name to Excel's 31-character limit
          safe_sheet_name <- substr(sheet_name, 1, 31)

          addWorksheet(wb, safe_sheet_name)
          writeDataTable(wb, sheet = safe_sheet_name, x = pheno_df,
                        tableName = gsub("[^A-Za-z0-9]", "", safe_sheet_name),
                        withFilter = TRUE)

          # Apply text wrapping
          if (wrap_text) {
            addStyle(wb, sheet = safe_sheet_name, style = wrap_style,
                    rows = 1:(nrow(pheno_df) + 1),
                    cols = 1:ncol(pheno_df),
                    gridExpand = TRUE)
          }

          # Auto-fit widths
          if (auto_width) {
            setColWidths(wb, sheet = safe_sheet_name,
                        cols = 1:ncol(pheno_df),
                        widths = "auto")
          }

          # Freeze header row
          if (freeze_panes) {
            freezePane(wb, sheet = safe_sheet_name, firstActiveRow = 2)
          }
        }
      }
    }
  }

  # Save workbook
  message("Creating Excel workbook: ", output_file)
  saveWorkbook(wb, output_file, overwrite = TRUE)

  invisible(output_file)
}


# =============================================================================
# Specialized Excel Functions (Study-Specific)
# =============================================================================

#' Create gene workbook with consolidated phenotype sheets
#'
#' Version of create_gene_workbook() that consolidates all phenotype/disease
#' data into two sheets: "AllMousePhenotypes" and "AllHumanDiseases" rather
#' than creating individual sheets per gene.
#'
#' @param gene_table Data frame with master gene table
#' @param output_file Character, path to output .xlsx file
#' @param mouse_phenotypes Optional data frame with mouse phenotype annotations
#' @param human_diseases Optional data frame with human disease associations
#' @param gene_id_col Character, gene ID column for joining (default "mouse_ensembl_id")
#' @param gene_symbol_col Character, gene symbol column (default "mouse_gene_symbol")
#' @param ... Additional arguments passed to create_gene_workbook()
#'
#' @return Invisibly returns the path to the created Excel file
#'
#' @details
#' This function creates consolidated phenotype sheets where each row represents
#' a gene-phenotype pair. This is more scalable than individual sheets per gene
#' and allows for easier filtering across all genes.
#'
#' @export
create_gene_workbook_consolidated <- function(gene_table,
                                               output_file,
                                               mouse_phenotypes = NULL,
                                               human_diseases = NULL,
                                               gene_id_col = "mouse_ensembl_id",
                                               gene_symbol_col = "mouse_gene_symbol",
                                               ...) {

  # Prepare phenotype tables list
  pheno_list <- list()

  if (!is.null(mouse_phenotypes) && nrow(mouse_phenotypes) > 0) {
    pheno_list[["AllMousePhenotypes"]] <- mouse_phenotypes
  }

  if (!is.null(human_diseases) && nrow(human_diseases) > 0) {
    pheno_list[["AllHumanDiseases"]] <- human_diseases
  }

  # Call main function
  create_gene_workbook(
    gene_table = gene_table,
    output_file = output_file,
    phenotype_tables = if (length(pheno_list) > 0) pheno_list else NULL,
    gene_id_col = gene_id_col,
    ...
  )
}


# =============================================================================
# Helper Functions for Preparing Phenotype Tables
# =============================================================================

#' Prepare mouse phenotype table from MGI/MouseMine data
#'
#' Reformats raw MouseMine phenotype data into clean table for Excel export.
#'
#' @param phenotype_data Data frame with MouseMine phenotype annotations
#' @param gene_col Character, name of gene column
#' @param term_id_col Character, name of ontology term ID column
#' @param term_name_col Character, name of ontology term name column
#' @param pubmed_col Character, name of PubMed ID column
#' @param description_col Character, name of description column (optional)
#'
#' @return Data frame formatted for Excel export
#'
#' @export
prepare_mouse_phenotype_table <- function(phenotype_data,
                                          gene_col = "mouse_gene_symbol",
                                          term_id_col = "OntologyAnnotation.ontologyTerm.identifier",
                                          term_name_col = "OntologyAnnotation.ontologyTerm.name",
                                          pubmed_col = "OntologyAnnotation.evidence.publications.pubMedId",
                                          description_col = NULL) {

  setDT(phenotype_data)

  # Select and rename columns
  selected_cols <- c(gene_col, term_id_col, term_name_col, pubmed_col)
  if (!is.null(description_col) && description_col %in% names(phenotype_data)) {
    selected_cols <- c(selected_cols, description_col)
  }

  pheno_table <- phenotype_data[, ..selected_cols]

  # Rename to cleaner names
  setnames(pheno_table,
           old = c(gene_col, term_id_col, term_name_col, pubmed_col),
           new = c("Gene", "OntologyTermID", "OntologyTermName", "PubMedID"),
           skip_absent = TRUE)

  if (!is.null(description_col) && description_col %in% names(phenotype_data)) {
    setnames(pheno_table, description_col, "Description", skip_absent = TRUE)
  }

  return(as.data.frame(pheno_table))
}


#' Prepare human disease association table from Open Targets
#'
#' Reformats raw Open Targets disease association data into clean table for Excel export.
#'
#' @param disease_data Data frame with Open Targets disease associations
#' @param gene_col Character, name of gene column (Ensembl ID or symbol)
#' @param disease_id_col Character, name of disease ID column
#' @param disease_name_col Character, name of disease name column
#' @param score_col Character, name of association score column
#'
#' @return Data frame formatted for Excel export
#'
#' @export
prepare_human_disease_table <- function(disease_data,
                                        gene_col = "human_ensembl_id",
                                        disease_id_col = "disease_id",
                                        disease_name_col = "disease_name",
                                        score_col = "association_score") {

  setDT(disease_data)

  selected_cols <- c(gene_col, disease_id_col, disease_name_col, score_col)
  disease_table <- disease_data[, ..selected_cols]

  # Rename to cleaner names
  setnames(disease_table,
           old = c(gene_col, disease_id_col, disease_name_col, score_col),
           new = c("Gene", "DiseaseID", "DiseaseName", "AssociationScore"),
           skip_absent = TRUE)

  # Round score to 3 decimal places
  if ("AssociationScore" %in% names(disease_table)) {
    disease_table[, AssociationScore := round(AssociationScore, 3)]
  }

  # Order by gene, then by score (descending)
  setorder(disease_table, Gene, -AssociationScore)

  return(as.data.frame(disease_table))
}


# =============================================================================
# Variant Table Generation
# =============================================================================

#' Create variant table for a locus
#'
#' Generates a CSV table of coding variants (SNPs/indels) within a locus.
#' Useful for including in locus packets.
#'
#' @param variants_in_locus Data frame with variant data
#' @param output_file Character, path to output .csv file
#' @param variant_id_col Character, name of variant ID column
#' @param chr_col Character, name of chromosome column
#' @param pos_col Character, name of position column
#' @param ref_col Character, name of reference allele column
#' @param alt_col Character, name of alternate allele column
#' @param consequence_col Character, name of consequence/effect column (optional)
#' @param gene_col Character, name of gene column (optional)
#'
#' @return Invisibly returns the path to the created CSV file
#'
#' @export
create_variant_table <- function(variants_in_locus,
                                  output_file,
                                  variant_id_col = "variant_id",
                                  chr_col = "chr",
                                  pos_col = "pos",
                                  ref_col = "ref",
                                  alt_col = "alt",
                                  consequence_col = NULL,
                                  gene_col = NULL) {

  if (nrow(variants_in_locus) == 0) {
    message("No variants found in locus. Skipping variant table.")
    return(invisible(NULL))
  }

  setDT(variants_in_locus)

  # Select columns
  selected_cols <- c(variant_id_col, chr_col, pos_col, ref_col, alt_col)
  if (!is.null(consequence_col) && consequence_col %in% names(variants_in_locus)) {
    selected_cols <- c(selected_cols, consequence_col)
  }
  if (!is.null(gene_col) && gene_col %in% names(variants_in_locus)) {
    selected_cols <- c(selected_cols, gene_col)
  }

  variant_table <- variants_in_locus[, ..selected_cols]

  # Write to CSV
  message("Creating variant table: ", output_file)
  fwrite(variant_table, output_file)

  invisible(output_file)
}
