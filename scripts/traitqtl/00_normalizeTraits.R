# --- Load Libraries, data, and args ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(optparse)
})

# Define the command-line arguments
option_list <- list(
  make_option(c("--input_raw"),  type = "character", help = "Path to input phenotype CSV file", metavar = "FILE"),
  make_option(c("--output_boxcox"), type = "character", help = "Path to output processed phenotype CSV file", metavar = "FILE"),
  make_option(c("--drug"),   type = "character", help = "Drug treatment to filter (e.g., Ctrl, Iso)", default = NULL)
)

# Parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$options$input_raw) || is.null(opt$options$output_boxcox)) {
  stop("Both --input_raw and --output_boxcox must be provided.")
}
if (!file.exists(opt$options$input_raw)) {
  stop("Input file does not exist: ", opt$options$input_raw)
}

# Read the input phenotype file
phenotypes <- read.csv(opt$options$input_raw)

# Filter the data based on the drug treatment (if requested and column exists)
if (!is.null(opt$options$drug)) {
  if (!"Drug" %in% names(phenotypes)) {
    stop("Column 'Drug' not found, but --drug was provided.")
  }
  phenotypes <- phenotypes |> dplyr::filter(Drug == opt$options$drug)
}

# List columns with measured or extrapolated phenotypes
traits <- phenotypes |> 
  dplyr::select(BW.day.0:CSA) |> 
  colnames()

# --- Data Preprocessing ---
# - Drop optional columns if present
# - Coerce trait columns to numeric 
# - Convert remaining character columns to factors
pheno_processed <- phenotypes |>
  dplyr::select(-any_of(c("Start.date", "End.date", "Notes"))) |>
  mutate(across(all_of(traits), ~ suppressWarnings(as.numeric(as.character(.))))) |>
  mutate(across(where(is.character), as.factor))

# --- Normalization (Box–Cox per trait per treatment) ---
# Notes:
# - Shifts traits with non-positive minima by +(-min + 1e-6)
# - Skips traits with all NA or zero variance
# - Uses MASS::boxcox(lm(y ~ 1), plotit = FALSE) to pick lambda
# - Applies log when lambda ~ 0

epsilon <- 1e-6
transform_log <- list()

for (trait in traits) {
  x <- pheno_processed[[trait]]

  # Skip if all NA
  if (all(is.na(x))) {
    warning(sprintf("Trait %s has all NA. Skipping.", trait))
    next
  }

  # Ensure numeric (already attempted above, but double-guard)
  if (!is.numeric(x)) {
    warning(sprintf("Trait %s is not numeric after coercion. Skipping.", trait))
    next
  }

  # Shift if min <= 0
  min_val <- suppressWarnings(min(x, na.rm = TRUE))
  shift_applied <- 0
  if (is.finite(min_val) && min_val <= 0) {
    shift_applied <- -min_val + epsilon
    x <- x + shift_applied
    message(sprintf("Trait %s had min=%.6g. Applied shift of %.6g.", trait, min_val, shift_applied))
  } else if (!is.finite(min_val)) {
    warning(sprintf("Trait %s minimum is not finite (all NA?). Skipping.", trait))
    next
  }

  # Non-NA vector after potential shift
  nn <- x[is.finite(x)]
  if (length(nn) < 2L) {
    warning(sprintf("Trait %s has <2 finite values. Skipping.", trait))
    next
  }
  if (sd(nn) == 0) {
    warning(sprintf("Trait %s has zero variance. Skipping Box–Cox.", trait))
    pheno_processed[[trait]] <- x
    transform_log[[trait]] <- list(shift = shift_applied, lambda = NA_real_, note = "zero variance")
    next
  }

  # Box–Cox lambda (no plotting). MASS::boxcox does NOT have 'quietly' arg.
  bc <- MASS::boxcox(lm(nn ~ 1), plotit = FALSE)
  lambda <- bc$x[which.max(bc$y)]

  # Apply transform to full vector x (already shifted if needed)
  if (abs(lambda) < 1e-6) {
    x_trans <- log(x)
  } else {
    x_trans <- (x^lambda - 1) / lambda
  }

  pheno_processed[[trait]] <- x_trans
  transform_log[[trait]] <- list(shift = shift_applied, lambda = lambda, note = "ok")
}

# Add a temp id for GWAS
pheno_processed <- pheno_processed %>%
  mutate(gwas_temp_id = paste0(Strain, Sex, Drug)) 

# --- Save Output ---

# Extract output file path from command-line arguments
output_file <- opt$options$output_boxcox

# Ensure output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  print(paste("Creating output directory:", output_dir))
  dir.create(output_dir, recursive = TRUE)
}

print(paste("Saving processed data to:", output_file))
write.csv(pheno_processed, file = output_file, row.names = F)