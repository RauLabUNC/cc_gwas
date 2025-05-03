# =============================================================================
# Harmonize QTL and eQTL Results from miQTL and pyLMM Analyses
#
# Description:
# This script reads output files from miQTL (interval-based) and pyLMM
# (SNP-based) analyses for both trait QTLs and eQTLs. It processes
# and formats these results into two harmonized tables: one for trait loci
# and one for eQTLs, saving them as CSV files.
#
# Author: Gemini
# Date: 2025-05-03
# =============================================================================

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
library(tidyr)
library(data.table) # Using data.table for efficient grouping of pyLMM SNPs

# --- 2. Define Parameters ---
# Significance thresholds
pyLMM_trait_pval_thresh <- 1e-5  # P-value threshold for pyLMM trait QTL SNPs
miQTL_eqtl_lod_thresh   <- 5.0   # LOD threshold for miQTL eQTLs
pyLMM_eqtl_pval_thresh  <- 1e-5  # P-value threshold for pyLMM eQTLs

# Grouping parameters for pyLMM trait QTL SNPs
pyLMM_grouping_dist <- 500000 # Max distance (bp) to group consecutive significant SNPs

# --- 3. Define File Paths ---
input_dir <- "data/processed/joinLoci"
output_dir <- "data/processed/joinLoci"

# Construct full file paths
miqtl_trait_loci_file <- file.path(input_dir, "all_significant_regions_summary.csv")
pylmm_trait_pvals_file <- file.path(input_dir, "Ctrl_pvals.csv")
miqtl_eqtl_file <- file.path(input_dir, "miQTL_output.csv")
pylmm_eqtl_file <- file.path(input_dir, "PyLMM_example_output.csv")

harmonized_trait_loci_outfile <- file.path(output_dir, "harmonized_trait_loci.csv")
harmonized_eqtl_outfile <- file.path(output_dir, "harmonized_eqtl.csv")

# --- 4. Load Data ---
miqtl_trait_loci_raw <- read.csv(miqtl_trait_loci_file)
pylmm_trait_pvals_raw <- read.csv(pylmm_trait_pvals_file)
miqtl_eqtl_raw <- read.csv(miqtl_eqtl_file, row.names = 1)
pylmm_eqtl_raw <- read.csv(pylmm_eqtl_file, row.names = 1)

# --- 5. Process miQTL Trait Loci ---
miqtl_trait_loci_processed <- miqtl_trait_loci_raw %>%
  dplyr::select(
      Trait = trait,
      Condition = drug,
      Chr = chr,
      Locus_Start_bp = lower_pos_lod_drop, # Will convert Mb -> bp later
      Peak_BP = peak_pos,                 # Will convert Mb -> bp later
      Locus_End_bp = upper_pos_lod_drop,   # Will convert Mb -> bp later
      Peak_LOD = max_lod,
      Normalization = norm # Keep normalization info if needed
    ) %>%
    mutate(
      Method = "miQTL-haplotype",
      # Convert positions from Mb (assumed) to bp
      Locus_Start_bp = as.numeric(Locus_Start_bp) * 1e6,
      Peak_BP = as.numeric(Peak_BP) * 1e6,
      Locus_End_bp = as.numeric(Locus_End_bp) * 1e6,
      Chr = as.character(Chr) # Ensure consistent Chr format
    ) %>%
    # Ensure columns are in a standard order
    dplyr::select(Method, Trait, Condition, Chr, Locus_Start_bp, Peak_BP, Locus_End_bp, Peak_LOD, Normalization)


# --- 6. Process pyLMM Trait Loci ---
pylmm_trait_loci_processed <- tibble() # Initialize empty tibble

# Check for expected SNP info columns
snp_info_cols <- c("Name", "Chr", "Pos_mm39") # Using mm39 position
if (!all(snp_info_cols %in% names(pylmm_trait_pvals_raw))) {
  stop(paste("pyLMM trait file missing expected SNP info columns:", paste(snp_info_cols[!snp_info_cols %in% names(pylmm_trait_pvals_raw)], collapse=", ")))
}

# Identify trait columns (assuming they are the ones remaining after SNP info)
trait_cols <- setdiff(names(pylmm_trait_pvals_raw),
                                c("Name", "Chr", "Pos_mm10", "Pos_mm39", "Main", "Alt", "rsID", "MAF"))
# Reshape from wide to long format
pylmm_trait_long <- pylmm_trait_pvals_raw %>%
  dplyr::select(all_of(snp_info_cols), all_of(potential_trait_cols)) %>%
  pivot_longer(
    cols = all_of(potential_trait_cols),
    names_to = "Trait",
    values_to = "P_VALUE"
  ) %>%
  filter(!is.na(P_VALUE)) %>% # Remove rows where P_VALUE might be NA after pivoting
  dplyr::rename(POS = Pos_mm39, # Use mm39 position
         CHR = Chr) %>%
  mutate(CHR = as.character(CHR))

# Filter significant SNPs
pylmm_sig_snps <- pylmm_trait_long %>%
  filter(P_VALUE < pyLMM_trait_pval_thresh)

if(nrow(pylmm_sig_snps) > 0) {
  # Use data.table for efficient grouping by distance
  setDT(pylmm_sig_snps)
  setorder(pylmm_sig_snps, Trait, CHR, POS) # Sort for grouping
  
  # Calculate difference in position with the previous SNP within each Trait/CHR group
  pylmm_sig_snps[, pos_diff := POS - shift(POS, type = "lag"), by = .(Trait, CHR)]
  
  # Assign a group ID based on the distance threshold
  pylmm_sig_snps[, group_id := cumsum(ifelse(is.na(pos_diff) | pos_diff > pyLMM_grouping_dist, 1, 0)), by = .(Trait, CHR)]
  
  # Summarize each group to define a locus
  pylmm_trait_loci_processed <- pylmm_sig_snps[, .(
    Locus_Start_bp = min(POS),
    Locus_End_bp = max(POS),
    Peak_BP = POS[which.min(P_VALUE)], # Position of the SNP with the minimum P-value
    Peak_Pval = min(P_VALUE)           # Minimum P-value in the group
  ), by = .(Trait, CHR, group_id)] %>%
    mutate(
      Peak_LOD = -log10(Peak_Pval),
      Method = "pyLMM-SNP",
      Condition = "Ctrl", # Assumed based on filename 'Ctrl_pvals.csv'
      Chr = as.character(CHR)
    ) %>%
    dplyr::select(Method, Trait, Condition, Chr, Locus_Start_bp, Peak_BP, Locus_End_bp, Peak_LOD) %>%
    as_tibble() # Convert back to tibble
  
  print(paste("Processed", nrow(pylmm_trait_loci_processed), "pyLMM trait loci from", nrow(pylmm_sig_snps), "significant SNP-trait associations."))
  
} else {
  print("No significant SNPs found in pyLMM trait file after filtering.")
}


# --- 7. Combine Trait Loci ---
# Use bind_rows, which handles differing columns by filling with NA
harmonized_trait_loci <- bind_rows(miqtl_trait_loci_processed, pylmm_trait_loci_processed)

print(paste("Total harmonized trait loci:", nrow(harmonized_trait_loci)))

# --- 8. Process miQTL eQTL ---
# Check for expected columns based on provided str() output
expected_cols <- c("trait", "lead.snp", "chr_region", "lead.snp.LOD", "treatment")
if (!all(expected_cols %in% names(miqtl_eqtl_raw))) {
  stop(paste("miQTL eQTL file missing expected columns:", paste(expected_cols[!expected_cols %in% names(miqtl_eqtl_raw)], collapse=", ")))
}

miqtl_eqtl_processed <- miqtl_eqtl_raw %>%
  filter(lead.snp.LOD >= miQTL_eqtl_lod_thresh) %>%
  # Extract Chr from chr_region (e.g., "1:3000355-197158220")
  mutate(Chr = str_extract(chr_region, "^[^:]+")) %>%
  dplyr::select(
    Gene_ID = trait,       # Gene ID is in 'trait' column for eQTL
    SNP_ID = lead.snp,     # Peak SNP identifier
    Chr,                   # Extracted chromosome
    Peak_LOD = lead.snp.LOD, # LOD score for the peak SNP
    Condition = treatment  # Condition/treatment group
  ) %>%
  mutate(
    Method = "miQTL-haplotype",
    # NOTE: Peak_BP (position of lead.snp) is not available in the provided str() columns.
    # Setting to NA. This needs to be obtained from another source if required.
    Peak_BP = NA_real_,
    Chr = as.character(Chr)
  ) %>%
  # Ensure standard column order (Peak_BP added, even if NA)
  dplyr::select(Method, Condition, Gene_ID, SNP_ID, Chr, Peak_BP, Peak_LOD)

print(paste("Processed", nrow(miqtl_eqtl_processed), "significant miQTL eQTLs."))


# --- 9. Process pyLMM eQTL ---
pylmm_eqtl_processed <- tibble() # Initialize

# Check for expected SNP info columns
snp_info_cols_eqtl <- c("Name", "Chr", "Pos_mm39") # Using mm39 position
if (!all(snp_info_cols_eqtl %in% names(pylmm_eqtl_raw))) {
  stop(paste("pyLMM eQTL file missing expected SNP info columns:", paste(snp_info_cols_eqtl[!snp_info_cols_eqtl %in% names(pylmm_eqtl_raw)], collapse=", ")))
}

# Identify gene columns (assuming they are the ones remaining after SNP info)
potential_gene_cols <- setdiff(names(pylmm_eqtl_raw),
                               c("Name", "Chr", "Pos_mm10", "Pos_mm39", "Main", "Alt", "rsID", "MAF"))
print(paste("Identified potential pyLMM gene columns:", paste(potential_gene_cols, collapse=", ")))

# Reshape from wide to long format
pylmm_eqtl_long <- pylmm_eqtl_raw %>%
  dplyr::select(all_of(snp_info_cols_eqtl), all_of(potential_gene_cols)) %>%
  pivot_longer(
    cols = all_of(potential_gene_cols),
    names_to = "Gene_ID", # Gene symbols are the column names
    values_to = "P_VALUE"
  ) %>%
  filter(!is.na(P_VALUE)) %>%
  dplyr::rename(POS = Pos_mm39, # Use mm39 position
         CHR = Chr,
         SNP_ID = Name) %>% # Use 'Name' column as SNP identifier
  mutate(CHR = as.character(CHR))

# Filter significant SNP-gene pairs
pylmm_eqtl_processed <- pylmm_eqtl_long %>%
  filter(P_VALUE <= pyLMM_eqtl_pval_thresh) %>%
  dplyr::select(
    SNP_ID,
    Chr = CHR,
    Peak_BP = POS,
    Gene_ID,
    Peak_Pval = P_VALUE
  ) %>%
  mutate(
    Peak_LOD = -log10(Peak_Pval),
    Method = "pyLMM-SNP",
    Condition = "Unknown", # Add condition if available, otherwise placeholder
    Chr = as.character(Chr),
    Peak_BP = as.numeric(Peak_BP)
  ) %>%
  # Ensure standard column order
  dplyr::select(Method, Condition, Gene_ID, SNP_ID, Chr, Peak_BP, Peak_LOD)

print(paste("Processed", nrow(pylmm_eqtl_processed), "significant pyLMM eQTL associations."))


# --- 10. Combine eQTL ---
harmonized_eqtl <- bind_rows(miqtl_eqtl_processed, pylmm_eqtl_processed)

print(paste("Total harmonized eQTLs:", nrow(harmonized_eqtl)))

# --- 11. Save Outputs ---
if (nrow(harmonized_trait_loci) > 0) {
  tryCatch({
    write_csv(harmonized_trait_loci, harmonized_trait_loci_outfile)
    print(paste("Harmonized trait loci saved to:", harmonized_trait_loci_outfile))
  }, error = function(e) {
    # Stop if saving fails, as it indicates a problem
    stop(paste("Error saving harmonized trait loci:", e$message))
  })
} else {
  print("No harmonized trait loci data to save.")
}

if (nrow(harmonized_eqtl) > 0) {
  tryCatch({
    write_csv(harmonized_eqtl, harmonized_eqtl_outfile)
    print(paste("Harmonized eQTLs saved to:", dplyr::))
  }, error = function(e) {
    # Stop if saving fails
    stop(paste("Error saving harmonized eQTLs:", e$message))
  })
} else {
  print("No harmonized eQTL data to save.")
}

print("Script finished.")
