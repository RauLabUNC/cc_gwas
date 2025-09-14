#this is mostly just a modified version of what is here https://github.com/gkeele/miqtl

#take a look at it for more explanations of what each step is doing
#remotes::install_github("gkeele/miqtl")
library(miqtl) 
library(tidyverse)
source("/proj/raulab/projects/cc_gwas/scripts/extras/scan_h2lmm.R")
genomecache <- "/proj/raulab/projects/cc_gwas/data/raw/genomes/Corrected_Haplotypes_Cache_CC_83024"

# Read in the first trailing argument for phenotype and treatment
args <- commandArgs(trailingOnly = TRUE)


#Load in phenotype data = gene expression ##BE VERY CAREFUL WITH WHICH COUNT FILE TO USE##
phenotypes <- read.csv("Data/Processed/ExpressionData/5d_sepVST_Info_250429.csv", check.names=F, row.names = 1) 
phenotypes$Drug <- ifelse(phenotypes$Drug == "Ctrl", 0, 1)
phenotypes$Sex <- ifelse(phenotypes$Sex == "M", 0, 1)



#Run eQTL on ctrl only, iso only, or both
if(args[2] == "control"){
  phenotypes <- phenotypes %>% filter(Drug == 0)
}else{if(args[2] == "iso"){
  phenotypes <- phenotypes %>% filter(Drug == 1)
  }else{phenotypes <- phenotypes}
}

# Run genome scan
  miqtl.rop.scan.scaled <- scan.h2lmm.test(
    genomecache = genomecache,
    data = phenotypes,
    pheno.id="SampleID",
    geno.id="Strain",
    formula = get(args[1]) ~ 1 + Sex,  #For this run, changed 0 to 1 in the linear model
    use.multi.impute = F,
    return.allele.effects = T,
    use.fix.par = T,
    chr = "all")


# Ensure the directory exists
output_dir <- file.path("Data/Scans",args[3], args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_file <- file.path(output_dir, paste0(args[1],"_scan_results.rds"))

# Save the output as an RDS file
saveRDS(miqtl.rop.scan.scaled, output_file)