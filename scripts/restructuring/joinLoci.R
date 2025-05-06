# =============================================================================
# Harmonize QTL and eQTL Results from miQTL and pyLMM Analyses
#
# Description:
# This script reads output files from miQTL (interval-based) and pyLMM
# (SNP-based) analyses for both trait QTLs and eQTLs. It processes
# and formats these results into two harmonized tables: one for trait loci
# and one for eQTLs, saving them as CSV files.

# =============================================================================

# --- 1. Load Libraries ---
library(dplyr)
library(readr)
library(tidyr)
library(data.table) # Using data.table for efficient grouping of pyLMM SNPs


#### Trait QTLs ####

## miQTL 
input_dir <- "data/processed/joinLoci/trait_qtl"
miqtl_trait_loci_file <- file.path(input_dir, "miQTL","all_significant_regions_summary.csv")
miqtl_trait_loci_raw <- read.csv(miqtl_trait_loci_file)


## Make reform to unify
miqlt_trait_loci_formatted <- miqtl_trait_loci_raw |> 
  mutate(Method = "miQTL",
         Analysis_Type = "TraitQTL",
         Dependent_Variable = trait,
         Treatment = drug,
         Chr = chr,
         Locus_Start_bp = upper_pos_lod_drop *10^6,
         Locus_End_bp = lower_pos_lod_drop*10^6,
         Peak_SNP_ID = Peak_SNP_ID,
         Peak_SNP_pos_bp = peak_pos*10^6,
         Lead_Strain = lead_strain,
         Peak_Significance_Value = max_lod,
         Significance_Metric = "LOD",
         Significance_Threshold = NA,
         Locus_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp,"_", 
                           Analysis_Type, "_",Dependent_Variable, "_", Treatment),
         Position_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp)) |>
  dplyr::select(Locus_ID, Position_ID, Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Lead_Strain, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)


## pyLMM 
pyLMM_trait_iso <- file.path(input_dir, "PyLMM","Iso_pvals.csv") |> read.csv()
pyLMM_trait_ctrl <- file.path(input_dir, "PyLMM","Ctrl_pvals.csv") |> read.csv()

trait_cols <- pyLMM_trait_iso |> 
  dplyr::select(BW.day.0:CSA) |> 
  colnames()

## Make reform to unify
pyLMM_trait_iso_formatted <- pyLMM_trait_iso |> 
  pivot_longer(
    cols = all_of(trait_cols),
    names_to = "Dependent_Variable",
    values_to = "Peak_Significance_Value"
  ) |> 
  mutate(Method = "PyLMM",
         Analysis_Type = "TraitQTL",
         Dependent_Variable = Dependent_Variable,
         Chr = Chr,
         Treatment = "Iso",
         Locus_Start_bp = Pos_mm39 - 10^6,
         Locus_End_bp = Pos_mm39 + 10^6,
         Peak_SNP_ID = Name,
         Peak_SNP_pos_bp = Pos_mm39,
         Significance_Metric = "P-value",
         Significance_Threshold = NA) |> 
  dplyr::select(Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)

pyLMM_trait_ctrl_formatted <- pyLMM_trait_ctrl |> 
  pivot_longer(
    cols = all_of(trait_cols),
    names_to = "Dependent_Variable",
    values_to = "Peak_Significance_Value"
  ) |> 
  mutate(Method = "PyLMM",
         Analysis_Type = "TraitQTL",
         Dependent_Variable = Dependent_Variable,
         Chr = Chr,
         Treatment = "Ctrl",
         Locus_Start_bp = Pos_mm39 - 5*10^5,
         Locus_End_bp = Pos_mm39 + 5*10^5,
         Peak_SNP_ID = Name,
         Peak_SNP_pos_bp = Pos_mm39,
         Significance_Metric = "P-value",
         Significance_Threshold = NA) |> 
  dplyr::select(Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)


# Bind the two tables
pyLMM_trait_formatted <- rbind(pyLMM_trait_ctrl_formatted, pyLMM_trait_iso_formatted) |> 
  mutate(Locus_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp,"_", 
                      Analysis_Type, "_",Dependent_Variable, "_", Treatment),
         Position_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp),
         Lead_Strain = NA
         
  ) |> 
  dplyr::select(Locus_ID, Position_ID, Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Lead_Strain, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)


#### eQTL data ####
input_dir <- "data/processed/joinLoci/eqtl"
miqtl_eqtl <- file.path(input_dir,"miQTL", "miQTL_output.csv") |> read.csv(row.names = 1)

## Make reform to unify
miqtl_eqtl_loci_formatted <- miqtl_eqtl |> 
  mutate(Method = "miQTL",
         Analysis_Type = "eQTL",
         Dependent_Variable = trait,
         Treatment = treatment,
         Chr = stringr::str_extract(chr_region, "^[^:]+"),
         Locus_Start_bp = as.numeric(stringr::str_extract(chr_region, "(?<=:)[0-9]+")),
         Locus_End_bp = as.numeric(stringr::str_extract(chr_region, "[0-9]+$")),
         Peak_SNP_ID = lead.snp,
         Peak_SNP_pos_bp = NA,
         Lead_Strain = lead.strain,
         Peak_Significance_Value = lead.snp.LOD,
         Significance_Metric = "LOD",
         Significance_Threshold = NA,
         Locus_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp,"_", 
                           Method, "_ ",Analysis_Type, "_",Dependent_Variable, "_", Treatment),
         Position_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp)) |>
  dplyr::select(Locus_ID, Position_ID, Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Lead_Strain, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)



## Make reform to unify
pyLMM_eqtl <- file.path(input_dir,"PyLMM", "PyLMM_example_output.csv") |> read.csv(row.names = 1)

not_traits <- pyLMM_eqtl |> 
  dplyr::select(Name:MAF) |> colnames()
traits <- colnames(pyLMM_eqtl)[!colnames(pyLMM_eqtl) %in% not_traits]

pyLMM_eqtl_formatted <- pyLMM_eqtl |> 
  pivot_longer(
    cols = all_of(traits),
    names_to = "Dependent_Variable",
    values_to = "Peak_Significance_Value"
  ) |> 
  mutate(Method = "PyLMM",
         Analysis_Type = "eQTL",
         Dependent_Variable = Dependent_Variable,
         Chr = factor(Chr),
         Treatment = factor("Ctrl"),
         Locus_Start_bp = Pos_mm39 - 10^6,
         Locus_End_bp = Pos_mm39 + 10^6,
         Peak_SNP_ID = Name,
         Peak_SNP_pos_bp = Pos_mm39,
         Lead_Strain = NA,
         Significance_Metric = "P-value",
         Significance_Threshold = NA,
         Locus_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp,"_", 
                           Method, "_ ",Analysis_Type, "_",Dependent_Variable, "_", Treatment),
         Position_ID = paste0(Chr, ":", Locus_Start_bp, "-",Locus_End_bp)) |> 
  dplyr::select(Locus_ID, Position_ID, Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Lead_Strain, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)



# --- Combine eQTL ---
trait_loci <- bind_rows(miqlt_trait_loci_formatted, pyLMM_trait_formatted)

exp_loci <- bind_rows(miqtl_eqtl_loci_formatted, pyLMM_eqtl_formatted)

## Save output
write.csv(trait_loci, "data/processed/joinLoci/relational_tables/traitLoci.csv", row.names = F)

## Save output
write.csv(exp_loci, "data/processed/joinLoci/relational_tables/expLoci.csv", row.names = F)
