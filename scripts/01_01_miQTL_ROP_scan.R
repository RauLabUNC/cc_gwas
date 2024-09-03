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
phenotypes <- read.csv("data/processed/phenotypes/mean_cc_panel_08_06_24.csv")

phenotypes_ctrl <- phenotypes |> 
  filter(!is.na(get(phenotype_of_interest)) & Drug_Clean == 0)

# Create a function to remove outliers for each given phenotype
remove_outliers <- function(data, variables) {
  for (var in variables) {
    # Calculate the lower and upper bounds for outliers
    q1 <- quantile(data[[var]], 0.25)
    q3 <- quantile(data[[var]], 0.75)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr
    
    # Remove outliers
    data <- data[!(data[[var]] < lower_bound | data[[var]] > upper_bound), ]
  }
  
  return(data)
}

phenotypes_ctrl <- remove_outliers(phenotypes_ctrl, phenotype_of_interest)


# Use phenotype_of_interest in the scan.h2lmm function

miqtl.rop.scan <- scan.h2lmm(
                            genomecache = genomecache,
                            data = phenotypes_ctrl,
                            pheno.id="Strain_Clean",
                            geno.id="Strain_Clean",
                            formula = get(phenotype_of_interest) ~ 1 + Sex_Clean,  #This is how you can incorporate co-variates
                            use.multi.impute = F,
                            return.allele.effects = T)

# Ensure the directory exists
output_dir <- "data/processed/scans"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_file <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_scan_results.rds"))

# Save the output as an RDS file
saveRDS(miqtl.rop.scan, output_file)
