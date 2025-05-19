#this is mostly just a modified version of what is here https://github.com/gkeele/miqtl

#take a look at it for more explanations of what each step is doing
#remotes::install_github("gkeele/miqtl")
library(miqtl) 
library(tidyverse)
source("scripts/extras/scan_h2lmm.R")
genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Read in the first trailing argument for phenotype and treatment
args <- commandArgs(trailingOnly = TRUE)

#args <- c("CSA", "iso")
if(args[2] == "control"){
  treatment <- "0"
}else{
  treatment <- "1"
}

#Load in phenotype data
phenotypes <- read.csv("data/processed/phenotypes/meanCenterScaledByTreat_03242025.csv")

phenotypes <- phenotypes |> filter(Drug_Binary == treatment)

# Run genome scan
miqtl.rop.scan.scaled <- scan.h2lmm.test(
  genomecache = genomecache,
  data = phenotypes,
  pheno.id="gwas_temp_id",
  geno.id="Strain",
  formula = get(args[1]) ~ 0 + Sex_Binary,  #This is how you can incorporate co-variates
  use.multi.impute = F,
  return.allele.effects = T,
  use.fix.par = T)

# Ensure the directory exists
output_dir <- file.path("data/processed/scans", args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_file <- file.path(output_dir, paste0(args[1], "_scan_results.rds"))

# Save the output as an RDS file
saveRDS(miqtl.rop.scan.scaled, output_file)
