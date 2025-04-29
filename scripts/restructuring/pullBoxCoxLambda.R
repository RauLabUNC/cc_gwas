# --- Load Libraries ---
library(MASS) # For boxcox()
library(tidyverse)


# --- Parameters ---
# Path to the original raw phenotype data CSV file
raw_pheno_file <- "data/raw/phenotypes/outliersTraits_03312025.csv" # VERIFY PATH

# Select the traits for which you want to calculate lambdas
# Using the full list from your previous run for thoroughness, adjust if needed
traits_to_check <- c(
  "BW.day.0", "BW.day.28", "TH", "LV", "RV", "LA", "RA", "Lung", "Liver", "Adrenal", 
  "THW.by.BW.0", "LVW.by.BW.0", "RVW.by.BW.0", "LAW.by.BW.0", "RAW.by.BW.0", 
  "LuW.by.BW.0", "LiW.by.BW.0", "AdrW.by.BW.0", "IVSd.0", "LVIDd.0", "LVPWd.0", 
  "IVSs.0", "LVIDs.0", "LVPWs.0", "HR.0", "EF.0", "FS.0", "LV.Mass.0", 
  "LV.Mass.Corrected.0", "LV.Vold.0", "LV.Vols.0", "LVIDd.by.LVIDs.0", 
  "LVIDd.by.LVIDs.28", "IVSd.28", "LVIDd.28", "LVPWd.28", "IVSs.28", "LVIDs.28", 
  "LVPWs.28", "HR.28", "EF.28", "FS.28", "LV.Mass.28", "LV.Mass.Corrected.28", 
  "LV.Vold.28", "LV.Vols.28", "delta.EF", "Wall.Thicknessd.0", 
  "Wall.Thicknessd.28", "Delta.Wall.Thickness.d", "Delta.LVIDd", "Delta.FS", 
  "Percent.Fibrosis", "CSA" 
) # CHANGE ME if needed

# --- Function to Calculate Optimal Lambda ---
# Simplified version of apply_boxcox, just returns lambda
get_optimal_lambda <- function(x, trait_name, group_name) {
  # Ensure numeric
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  
  # Handle non-positive values with a shift
  min_val <- min(x_num, na.rm = TRUE)
  shift <- 0 
  
  if (is.finite(min_val) && min_val <= 0) {
    shift <- -min_val + 1e-6
    # warning(paste("Trait", trait_name, "in group", group_name, "has non-positive min (", min_val, "). Applying shift:", shift))
    x_num <- x_num + shift
  } else if (!is.finite(min_val)) {
    warning(paste("Trait", trait_name, "in group", group_name, "min value not finite."))
    return(NA_real_) # Return NA if min is not finite
  }
  
  # Extract non-NA data
  non_na_data <- x_num[!is.na(x_num)]
  
  if (length(non_na_data) > 0 && all(non_na_data > 0, na.rm = TRUE)) {
    # Use tryCatch to handle potential errors in boxcox itself
    bc_result <- tryCatch({
      MASS::boxcox(non_na_data ~ 1, plotit = FALSE, quietly = TRUE)
    }, error = function(e) {
      warning(paste("Error during boxcox calculation for trait", trait_name, "in group", group_name, ":", e$message))
      return(NULL) # Return NULL if boxcox fails
    })
    
    if (!is.null(bc_result)) {
      lambda <- bc_result$x[which.max(bc_result$y)]
      return(lambda)
    } else {
      return(NA_real_) # Return NA if boxcox calculation failed
    }
    
  } else {
    warning(paste("Trait", trait_name, "in group", group_name, "has no positive data after shift or only NAs."))
    return(NA_real_) # Return NA if no valid data
  }
}

# --- Load Data ---
print(paste("Loading raw data from:", raw_pheno_file))
pheno_raw <- read.csv(raw_pheno_file)

# --- Prepare Data ---
# Ensure traits are numeric, filter for Ctrl/Iso
pheno_prep <- pheno_raw %>%
  select(SampleID, Drug, all_of(traits_to_check)) %>%
  mutate(across(all_of(traits_to_check), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  filter(Drug %in% c("Ctrl", "Iso"))

# --- Calculate Lambdas for each group ---
lambda_results <- list()

for (trait in traits_to_check) {
  print(paste("Calculating lambdas for trait:", trait))
  
  # Calculate for Ctrl group
  ctrl_data <- pheno_prep %>% filter(Drug == "Ctrl") %>% pull(!!sym(trait))
  lambda_ctrl <- get_optimal_lambda(ctrl_data, trait, "Ctrl")
  
  # Calculate for Iso group
  iso_data <- pheno_prep %>% filter(Drug == "Iso") %>% pull(!!sym(trait))
  lambda_iso <- get_optimal_lambda(iso_data, trait, "Iso")
  
  # Store results
  lambda_results[[trait]] <- tibble(
    Trait = trait,
    Lambda_Ctrl = lambda_ctrl,
    Lambda_Iso = lambda_iso
  )
}

# --- Combine and Display Results ---
final_lambdas <- bind_rows(lambda_results)

print("Optimal Lambdas Calculated Separately for Ctrl and Iso Groups:")
print(final_lambdas, n = nrow(final_lambdas)) # Print all rows

final_lambdas |> 
  mutate(delta = Lambda_Ctrl - Lambda_Iso) |> 
  arrange(Lambda_Ctrl)
# Optional: Save results to a file
# write.csv(final_lambdas, "results/split_lambdas_by_drug.csv", row.names = FALSE)

