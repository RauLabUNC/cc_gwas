#this is mostly just a modified version of what is here https://github.com/gkeele/miqtl

#take a look at it for more explanations of what each step is doing
#remotes::install_github("gkeele/miqtl")
library(miqtl) 
library(tidyverse)

genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]

#So, I've modified this file such that it has a column (strain_clean) that has
#just the strain names, nothing more.
phenotypes <- read.csv("data/processed/phenotypes/no_outliers_cc_panel_08_06_24.csv")

phenotypes_ctrl <- phenotypes |> 
  filter(!is.na(get(phenotype_of_interest)))

# Make co-variates categorical
phenotypes_ctrl <- phenotypes_ctrl |>
  mutate(Sex_Clean = as.factor(Sex_Clean),
         Drug_Clean = as.factor(Drug_Clean))

# Use phenotype_of_interest in the scan.h2lmm function
miqtl.rop.scan <- scan.h2lmm(
                            genomecache = genomecache,
                            data = phenotypes_ctrl,
                            pheno.id="Strain_Clean",
                            geno.id="Strain_Clean",
                            formula = get(phenotype_of_interest) ~ 1 + Sex_Clean + Drug_Clean,  #This is how you can incorporate co-variates
                            use.multi.impute = F,
                            return.allele.effects = T, 
                            use.fix.par = T)

# Ensure the directory exists
output_dir <- "data/processed/scans"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_file <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_scan_results.rds"))

# Save the output as an RDS file
saveRDS(miqtl.rop.scan, output_file)
