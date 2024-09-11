#Permutation.  This will get us the appropriate GWAS significance cutoff.  We 
# *probably* only have to do this once per condition we try - basically, as long as 
# the number of SNPs and Strains tested don't change and the nature of the data
# doesn't change much either (so ctrl vs ISO, for example), things should be fairly stable.
# but if this is all running in the background on the cluster and its fast enough... might as well
# do it every time :) 

# Load libraries
library(miqtl)
library(tidyverse)
# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]

# Load the data
scan_file <- file.path("data/processed/scans", paste0(as.character(phenotype_of_interest), "_scan_results.rds"))
scan <- readRDS(scan_file)

#So, I've modified this file such that it has a column (strain_clean) that has
#just the strain names, nothing more.
phenotypes <- read.csv("data/processed/phenotypes/no_outliers_cc_panel_08_06_24.csv")


phenotypes_ctrl <- phenotypes |> 
  mutate(pheno.id = paste(Strain_Clean, Sex_Clean, Drug_Clean, sep = "_")) |> 
  filter(!is.na(get(phenotype_of_interest)))

# Load genome cache
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Permute phenotypes and run scans on each
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan, 
                                                      method = "permutation", num.samples = 5,
                                                      use.BLUP = T, model.type = "null")

permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes_ctrl,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1)

permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, 
                                        percentile = 0.95) # 10.1186/s40168-023-01552-8 uses this percentile


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