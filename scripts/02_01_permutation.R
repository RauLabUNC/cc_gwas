# Load libraries
library(miqtl)
library(tidyverse)
# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]

# Load the data
scan_file <- file.path("data/processed/scans", paste0(as.character(phenotype_of_interest), "_scan_results.rds"))
scan <- readRDS(scan_file)

# Load phenotypes
phenotypes <- read.csv("data/processed/phenotypes/no_outliers_cc_panel_08_06_24.csv")

# Remove NAs, summarize groups to means, and scale
scaled_phenotypes <- phenotypes |> 
  filter(!is.na(get(phenotype_of_interest))) |>
  group_by(Strain_Clean, Drug_Clean, Sex_Clean, pheno.id) |> 
  summarize(across(BW.day.0:Percent.Fibrosis, ~mean(as.numeric(.), na.rm = TRUE))) |>
  ungroup() |> 
  mutate(across(BW.day.0:Percent.Fibrosis, ~ (. - mean(., na.rm = T))/sd(., na.rm = T))) |> 
  as.data.frame()


# Load genome cache
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Permute phenotypes and run scans on each
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan, 
                                                      method = "permutation", num.samples = 20,
                                                      use.BLUP = T, model.type = "null")

permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes_ctrl,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1)

permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, 
                                        percentile = 0.9) # 10.1186/s40168-023-01552-8 uses 85 percentile


# Ensure the directory exists
output_dir <- "data/processed/scan_thresholds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_threshold <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_threshold.rds"))
output_scan <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_scan.rds"))

# Save the threshold and scan as RDS files
saveRDS(permute_threshold, output_threshold)
saveRDS(permuted_scans, output_scan)