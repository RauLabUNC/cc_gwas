# --- Load Libraries, data, and args ---
library(tidyverse)
library(MASS) 
library(optparse)

# Define the command-line arguments
# Define the command-line arguments
option_list <- list(
  make_option(c("--input"), type = "character", help = "Path to input phenotype file"),
  make_option(c("--output"), type = "character", help = "Path to output processed phenotype file"),
  make_option(c("--normalization"), type = "character", help = "Normalization method (e.g., zscore, boxcox)"),
  make_option(c("--aggregation"), type = "character", help = "Aggregation method (e.g., individual, mean)"),
  make_option(c("--drug"), type = "character", help = "Drug treatment to filter (e.g., Ctrl, Iso)")
)
# Parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Print the parsed arguments (for debugging)
print(opt)

# Read the input phenotype file
phenotypes <- read.csv(opt$input)

# Filter the data based on the drug treatment
if (!is.null(opt$drug)) {
  phenotypes <- phenotypes |> filter(Drug == opt$drug)
}

# List columns with measured or extrapolated phenotypes
traits <- phenotypes |> 
  dplyr::select(BW.day.0:CSA) |> 
  colnames()

# --- Data Preprocessing ---

# Select columns and set types
pheno_processed <- phenotypes %>%
  dplyr::select(-Start.date, -End.date) |> 
  mutate(across(all_of(traits), ~ suppressWarnings(as.numeric(as.character(.))))) |> 
  # Convert character columns (like group_columns if they are character) to factors
  mutate(across(where(is.character), as.factor)) 

# --- Normalization ---

# Extract normalization method from command-line arguments
normalization_method <- opt$normalization

if (normalization_method == "zscore") {
  pheno_processed <- pheno_processed %>%
    mutate(across(all_of(traits), ~ scale(.)[, 1])) # scale() returns matrix, take first col
} else if (normalization_method == "boxcox") {
  # Loop through each specified trait for Box-Cox transformation
  for (trait in traits) {
    current_trait_data <- pheno_processed[[trait]]
    
    # Check for non-positive values using the overall minimum
    min_val <- min(current_trait_data, na.rm = TRUE)
    
    # Initialize shift_applied flag for this trait
    shift_applied_value <- 0 
    
    if (is.finite(min_val) && min_val <= 0) {
      # Calculate shift needed to make the minimum value slightly positive (e.g., 1e-6)
      shift <- -min_val + 1e-6 
      warning(paste("Trait", trait, "has non-positive minimum value (", min_val, "). Applying shift of", shift, "to ensure all values are positive."))
      pheno_processed[[trait]] <- current_trait_data + shift
      shift_applied_value <- shift # Record the shift applied
    } else if (!is.finite(min_val)) {
      warning(paste("Trait", trait, "minimum value is not finite (all NA?). Skipping shift."))
    }
    
    # Extract non-NA data for Box-Cox calculation *after potential shift*
    non_na_data <- pheno_processed[[trait]][!is.na(pheno_processed[[trait]])]
    
    # Check again if data is now positive before applying boxcox
  if (length(non_na_data) > 0 && all(non_na_data > 0, na.rm = TRUE)) { # Added na.rm=TRUE for safety
      # Perform Box-Cox calculation (requires a model, ~ 1 is simplest)
      # Suppress plot generation
      bc_result <- MASS::boxcox(lm(non_na_data ~ 1), plotit = FALSE, quietly = TRUE) 
      lambda <- bc_result$x[which.max(bc_result$y)] 
      
      # Apply the transformation
      if (abs(lambda) < 1e-6) { # Check if lambda is effectively zero (log transform)
        pheno_processed[[trait]] <- log(pheno_processed[[trait]])
      } else {
        pheno_processed[[trait]] <- (pheno_processed[[trait]]^lambda - 1) / lambda
      }
     }  
    }
   } 

# --- Aggregation ---
aggregation_method <- opt$aggregation

# --- Aggregation ---
if (aggregation_method == "mean") {
  pheno_processed <- pheno_processed %>%
    group_by(Drug, Sex, Strain) %>%
    # Calculate mean for numeric traits, handle potential NAs
    summarise(across(all_of(traits), ~ mean(., na.rm = TRUE)), 
              .groups = 'drop') }# Drop grouping structure after summarizing

# Add a temp id for GWAS
pheno_processed <- pheno_processed %>%
  mutate(gwas_temp_id = paste0(Strain, Sex, Drug)) 
# --- Save Output ---

# Extract output file path from command-line arguments
output_file <- opt$output

# Ensure output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  print(paste("Creating output directory:", output_dir))
  dir.create(output_dir, recursive = TRUE)
}

print(paste("Saving processed data to:", output_file))
write.csv(pheno_processed, file = output_file, row.names = F)
