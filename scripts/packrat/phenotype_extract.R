#' PackRat: Phenotype Keyword Extraction Functions
#'
#' Functions for extracting and filtering phenotype/disease associations
#' based on user-defined keyword patterns. Customizable for any phenotype class.
#'
#' @author Brian Gural, Anh Luu, Todd Kimball, Christoph Rau
#' @references Gural et al. (in preparation)

# =============================================================================
# Keyword-Based Phenotype Extraction
# =============================================================================

#' Extract phenotype terms matching keywords from a summary string
#'
#' @param phenotype_summary Character string with comma-separated phenotype terms
#' @param keywords Character vector of keywords/patterns to search for
#' @param ignore_case Logical, case-insensitive matching (default TRUE)
#' @param word_boundary Logical, require word boundaries for matches (default TRUE)
#' @param return_unique Logical, return unique matches only (default TRUE)
#'
#' @return Character string of matching phenotypes (comma-separated) or
#'   "No [type] phenotypes identified" if none found
#'
#' @details Uses regex pattern matching to identify phenotype terms containing
#'   any of the provided keywords. Designed for extracting disease-specific
#'   phenotypes (e.g., cardiac, neurological, metabolic) from comprehensive
#'   phenotype summaries.
#'
#' @examples
#' # Cardiac phenotype extraction
#' cardiac_keywords <- c("heart", "cardiac", "myocardium", "coronary")
#' pheno_string <- "abnormal heart morphology, kidney defect, cardiac hypertrophy"
#' extract_phenotypes(pheno_string, cardiac_keywords)
#' # Returns: "abnormal heart morphology, cardiac hypertrophy"
#'
#' # Neurological phenotype extraction
#' neuro_keywords <- c("brain", "neural", "cognitive", "seizure")
#' pheno_string <- "normal behavior, seizure susceptibility, heart defect"
#' extract_phenotypes(pheno_string, neuro_keywords)
#' # Returns: "seizure susceptibility"
#'
#' @export
extract_phenotypes <- function(phenotype_summary,
                               keywords,
                               ignore_case = TRUE,
                               word_boundary = TRUE,
                               return_unique = TRUE) {

  # Handle missing or empty input
  if (is.na(phenotype_summary) || phenotype_summary == "") {
    return("No phenotypes identified")
  }

  # Build regex pattern from keywords
  if (word_boundary) {
    pattern <- paste0("\\b(", paste(keywords, collapse = "|"), ")\\b")
  } else {
    pattern <- paste0("(", paste(keywords, collapse = "|"), ")")
  }

  # Split phenotype summary by delimiter (assumes comma-separated)
  all_phenotypes <- trimws(unlist(strsplit(phenotype_summary, ",")))

  # Filter for matching phenotypes
  matching_phenotypes <- all_phenotypes[grepl(pattern, all_phenotypes,
                                              ignore.case = ignore_case)]

  # Return results
  if (length(matching_phenotypes) == 0) {
    return("No phenotypes identified")
  } else {
    if (return_unique) {
      matching_phenotypes <- unique(matching_phenotypes)
    }
    return(paste(matching_phenotypes, collapse = ", "))
  }
}


#' Extract human disease associations matching keywords
#'
#' Wrapper for extract_phenotypes() specific to human disease data.
#'
#' @param disease_summary Character string with comma-separated disease terms
#' @param keywords Character vector of disease-related keywords
#' @param ... Additional arguments passed to extract_phenotypes()
#'
#' @return Character string of matching diseases or "No traits identified"
#'
#' @examples
#' cardiac_keywords <- c("heart", "cardiac", "coronary", "myocard")
#' disease_string <- "coronary artery disease, type 2 diabetes, hypertension"
#' extract_human_traits(disease_string, cardiac_keywords)
#'
#' @export
extract_human_traits <- function(disease_summary, keywords, ...) {
  result <- extract_phenotypes(disease_summary, keywords, ...)
  if (result == "No phenotypes identified") {
    return("No traits identified")
  }
  return(result)
}


#' Extract mouse phenotype associations matching keywords
#'
#' Wrapper for extract_phenotypes() specific to mouse phenotype data (e.g., MGI).
#'
#' @param phenotype_summary Character string with comma-separated phenotype terms
#' @param keywords Character vector of phenotype-related keywords
#' @param ... Additional arguments passed to extract_phenotypes()
#'
#' @return Character string of matching phenotypes or "No phenotypes identified"
#'
#' @examples
#' cardiac_keywords <- c("heart", "cardiac", "myocardium", "ventricular")
#' pheno_string <- "abnormal heart ventricle morphology, increased body weight"
#' extract_mouse_phenotypes(pheno_string, cardiac_keywords)
#'
#' @export
extract_mouse_phenotypes <- function(phenotype_summary, keywords, ...) {
  extract_phenotypes(phenotype_summary, keywords, ...)
}


# =============================================================================
# Formatted Output Functions (for reports/summaries)
# =============================================================================

#' Format gene-phenotype associations for human-only genes
#'
#' Creates formatted text output listing genes and their associated human traits.
#'
#' @param genes_df Data frame with columns 'Mouse Gene Symbol' and 'Human Disease Summary'
#' @param keywords Character vector of keywords to filter diseases
#' @param gene_col Name of column containing gene symbols (default "Mouse Gene Symbol")
#' @param pheno_col Name of column containing phenotype summary (default "Human Disease Summary")
#'
#' @return Character string with formatted gene:trait pairs (newline-separated)
#'
#' @examples
#' genes <- data.frame(
#'   `Mouse Gene Symbol` = c("Abcb10", "Acsl5"),
#'   `Human Disease Summary` = c("sepsis, cardiac arrest", "diabetes, obesity"),
#'   check.names = FALSE
#' )
#' cardiac_keywords <- c("cardiac", "heart")
#' format_human_only(genes, cardiac_keywords)
#' # Returns: "Abcb10: cardiac arrest"
#'
#' @export
format_human_only <- function(genes_df,
                               keywords,
                               gene_col = "Mouse Gene Symbol",
                               pheno_col = "Human Disease Summary") {

  if (nrow(genes_df) == 0) return("None identified")

  result <- character(nrow(genes_df))
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df[[gene_col]][i]
    traits <- extract_human_traits(genes_df[[pheno_col]][i], keywords)
    result[i] <- paste0(gene, ": ", traits)
  }
  return(paste(result, collapse = "\n        "))
}


#' Format gene-phenotype associations for mouse-only genes
#'
#' Creates formatted text output listing genes and their associated mouse phenotypes.
#'
#' @param genes_df Data frame with columns for gene symbols and phenotype summary
#' @param keywords Character vector of keywords to filter phenotypes
#' @param gene_col Name of column containing gene symbols
#' @param pheno_col Name of column containing phenotype summary
#'
#' @return Character string with formatted gene:phenotype pairs (newline-separated)
#'
#' @export
format_mouse_only <- function(genes_df,
                               keywords,
                               gene_col = "Mouse Gene Symbol",
                               pheno_col = "Mouse Phenotype Summary (MGI)") {

  if (nrow(genes_df) == 0) return("None identified")

  result <- character(nrow(genes_df))
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df[[gene_col]][i]
    phenotypes <- extract_mouse_phenotypes(genes_df[[pheno_col]][i], keywords)
    result[i] <- paste0(gene, ": ", phenotypes)
  }
  return(paste(result, collapse = "\n        "))
}


#' Format gene-phenotype associations for genes with both human and mouse evidence
#'
#' Creates formatted text output with hierarchical structure showing both human
#' traits and mouse phenotypes for each gene.
#'
#' @param genes_df Data frame with gene symbols, human disease, and mouse phenotype columns
#' @param keywords Character vector of keywords to filter phenotypes
#' @param gene_col Name of column containing gene symbols
#' @param human_col Name of column containing human disease summary
#' @param mouse_col Name of column containing mouse phenotype summary
#'
#' @return Character string with formatted gene:trait/phenotype pairs (newline-separated)
#'
#' @export
format_both <- function(genes_df,
                        keywords,
                        gene_col = "Mouse Gene Symbol",
                        human_col = "Human Disease Summary",
                        mouse_col = "Mouse Phenotype Summary (MGI)") {

  if (nrow(genes_df) == 0) return("None identified")

  result <- character(nrow(genes_df))
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df[[gene_col]][i]
    human_traits <- extract_human_traits(genes_df[[human_col]][i], keywords)
    mouse_phenotypes <- extract_mouse_phenotypes(genes_df[[mouse_col]][i], keywords)
    result[i] <- paste0(gene, ":\n            Human: ", human_traits,
                       "\n            Mouse: ", mouse_phenotypes)
  }
  return(paste(result, collapse = "\n\n        "))
}


# =============================================================================
# Predefined Keyword Sets (Examples)
# =============================================================================

#' Cardiac/cardiovascular phenotype keywords
#'
#' Comprehensive set of keywords for identifying cardiac-related phenotypes
#' in both human disease databases and mouse phenotype ontologies.
#'
#' @return Character vector of cardiac-related keywords
#' @export
cardiac_keywords <- function() {
  c(
    # Basic cardiac terms
    "heart", "cardi", "coronary", "atri", "ventric", "myocard",
    # Heart function issues
    "hypertroph", "fibril", "valv", "rhythm", "arrhythm", "tachycard", "bradycard",
    # Circulation and vascular
    "blood pressure", "pulse", "aort", "angina", "stroke", "hypertens", "thromb", "pulse pressure",
    # Ischemia and infarction
    "infarct", "ischem",
    # Mouse-specific structure terms
    "heart weight", "heart atrium", "heart ventricle", "heart morphology",
    "atrioventricular", "pericardial", "ventricular septal defect",
    "atrial", "thin myocardium", "common atrioventricular valve",
    # Vasculature terms
    "coronary vessel", "blood vessel", "vascular", "vasculature",
    "lymphatic vessel", "yolk sac vascular", "arch artery", "pharyngeal arch artery",
    # Circulation terms
    "heart failure", "circulat", "hemorrhage", "congest", "blood circulation",
    # Heart function terms
    "response of heart", "cardiac muscle", "myocardium layer", "myocardial fiber",
    # Abbreviations
    "LV\\."
  )
}


#' Example: Metabolic phenotype keywords
#'
#' Example keyword set for metabolic phenotypes. Users can define their own.
#'
#' @return Character vector of metabolic keywords
#' @export
metabolic_keywords <- function() {
  c("glucose", "insulin", "diabetes", "metabolic", "lipid", "cholesterol",
    "triglyceride", "obesity", "adipos", "hepatic steatosis", "NAFLD")
}


#' Example: Neurological phenotype keywords
#'
#' Example keyword set for neurological phenotypes. Users can define their own.
#'
#' @return Character vector of neurological keywords
#' @export
neurological_keywords <- function() {
  c("brain", "neural", "neuron", "cognitive", "behavior", "seizure",
    "epilep", "Alzheimer", "Parkinson", "dementia", "memory", "learning",
    "anxiety", "depression", "motor", "cerebral", "hippocampus")
}
