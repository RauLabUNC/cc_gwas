#this is mostly just a modified version of what is here https://github.com/gkeele/miqtl

#take a look at it for more explanations of what each step is doing
#remotes::install_github("gkeele/miqtl")
library(miqtl)
#"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"
genomecache <- "data/raw/genomes/CC_Genome_Cache_Clean_w_Founders"
genomecache2 <- "data/raw/genomes/CC_Genome_Cache"

# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
phenotype_of_interest <- args[1]
#This will simulate phenotypes.  It only works without founder lines.

#phenotypes <- sim.CC.data(genomecache = genomecache2, 
 #                         num.lines = 40, 
  #                        num.sim = 2, 
   #                       num.replicates = 2, 
    #                      qtl.effect.size = 0.8)


#So, I've modified this file such that it has a column (strain_clean) that has
#just the strain names, nothing more.

phenotypes <- read.csv("data/raw/phenotypes/full_cc_panel_data_04_16_24.csv")

phenotypes_ctrl <- phenotypes[phenotypes$Drug=="Ctrl",]

# Make column of interest into numeric type
phenotypes_ctrl[,phenotype_of_interest] <- as.numeric(phenotypes_ctrl[,phenotype_of_interest])

# QTL scan using ROP
# Function to return phenotype of interest without quotes

# Use phenotype_of_interest in the scan.h2lmm function
miqtl.rop.scan = scan.h2lmm(genomecache = genomecache,
                            data = phenotypes_ctrl,
                            pheno.id="Strain_Clean",
                            geno.id="Strain_Clean",
                            formula = get(phenotype_of_interest) ~ 1, #This is how you can incorporate covariates like sex or... other things.
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
