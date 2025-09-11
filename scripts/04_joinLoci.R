# =============================================================================
# Harmonize QTL and eQTL Results from miQTL Analyses
#
# Description:
# This script reads output files from miQTL analyses for both trait QTLs 
# and eQTLs. It processes and formats these results into two harmonized 
# tables: one for trait loci and one for eQTLs, saving them as CSV files.

# =============================================================================

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
library(tidyr)


# --- Trait QTLs ---

## miQTL trait loci
input_dir <- "data/processed/joinLoci/trait_qtl"
miqtl_trait_loci_file <- file.path(input_dir, "miQTL","all_significant_regions_summary.csv")
miqtl_trait_loci_raw <- read.csv(miqtl_trait_loci_file)

## Format trait QTL data
miqtl_trait_loci_formatted <- miqtl_trait_loci_raw |> 
  mutate(Method = "miQTL",
         Analysis_Type = "TraitQTL",
         Dependent_Variable = trait,
         Treatment = drug,
         Chr = chr,
         Locus_Start_bp = upper_pos_lod_drop * 10^6,
         Locus_End_bp = lower_pos_lod_drop * 10^6,
         Peak_SNP_ID = Peak_SNP_ID,
         Peak_SNP_pos_bp = peak_pos * 10^6,
         Lead_Strain = lead_strain,
         Peak_Significance_Value = max_lod,
         Significance_Metric = "LOD",
         Significance_Threshold = NA,
         Locus_ID = paste0(Chr, ":", Locus_Start_bp, "-", Locus_End_bp, "_", 
                           Analysis_Type, "_", Dependent_Variable, "_", Treatment),
         Position_ID = paste0(Chr, ":", Locus_Start_bp, "-", Locus_End_bp)) |>
  dplyr::select(Locus_ID, Position_ID, Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Lead_Strain, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)

# --- Save Outputs ---
trait_loci <- miqtl_trait_loci_formatted

if(!dir.exists("data/processed/joinLoci/relational_tables")){
  dir.create("data/processed/joinLoci/relational_tables", recursive = TRUE)
}

## Save output
write.csv(trait_loci, "data/processed/joinLoci/relational_tables/traitLoci.csv", row.names = F)
