#this is mostly just a modified version of what is here https://github.com/gkeele/miqtl

#take a look at it for more explanations of what each step is doing
#remotes::install_github("gkeele/miqtl")
library(miqtl) 
library(tidyverse)

genomecache <- "data/raw/genomes/haplotype_cache_cc_083024"

# Read in the first trailing argument for phenotype and treatment
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]

if(args[2] == "control"){
  treatment <- "0"
}else{
  treatment <- "1"
}
#Load in phenotype data
phenotypes <- read.csv("data/processed/phenotypes/no_outliers_cc_panel_08_06_24.csv")

# Remove NAs, summarize groups to means, and scale
scaled_phenotypes <- phenotypes |> 
  filter(!is.na(get(phenotype_of_interest)) & Drug_Clean == treatment) |>
  group_by(Strain_Clean, Drug_Clean, Sex_Clean, pheno.id) |> 
  summarize(across(BW.day.0:Percent.Fibrosis, ~mean(as.numeric(.), na.rm = TRUE))) |>
  ungroup() |> 
  mutate(across(BW.day.0:Percent.Fibrosis, ~ (. - mean(., na.rm = T))/sd(., na.rm = T))) |> 
  as.data.frame()

# Run genome scan
miqtl.rop.scan.scaled <- scan.h2lmm(
  genomecache = genomecache,
  data = scaled_phenotypes,
  pheno.id="pheno.id",
  geno.id="Strain_Clean",
  formula = get(phenotype_of_interest) ~ 0 + Sex_Clean,  #This is how you can incorporate co-variates
  use.multi.impute = F,
  return.allele.effects = T,
  use.fix.par = T)

# Ensure the directory exists
output_dir <- file.path("data/processed/scans", args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_file <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_scan_results.rds"))

# Save the output as an RDS file
saveRDS(miqtl.rop.scan.scaled, output_file)
