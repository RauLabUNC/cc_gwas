# Load libraries
library(miqtl)
library(tidyverse)
# Read in the first trailing argument for phenotype and treatment
args <- commandArgs(trailingOnly = TRUE)


########################################

# Load the data
scan_file <- file.path("Data/Scans", args[3], args[2], paste0(args[1],"_scan_results.rds"))
scan <- readRDS(scan_file)

# Load phenotypes ###WATCH WHICH FILE TO USE###
phenotypes <- read.csv("Data/Processed/ExpressionData/5d_sepVST_Info_250429.csv", check.names=F, row.names = 1)



# Filter for drug group(s) we want
if(args[2] == "control"){
  phenotypes <- phenotypes %>% filter(Drug == "Ctrl")
}else{if(args[2] == "iso"){
  phenotypes <- phenotypes %>% filter(Drug == "Iso")
}else{phenotypes <- phenotypes}
}


# Load genome cache
genomecache <- "/proj/raulab/projects/cc_gwas/data/raw/genomes/Corrected_Haplotypes_Cache_CC_83024"



#Get the gene chromosome number
gene_pos = read.csv("Data/Raw/genes_mouse.csv") %>%
  filter(mouse_gene_symbol == args[1])
gene_chromosome <- gene_pos$chr


# Permute phenotypes and run scans on each
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan,
                                                      subsample.chr = gene_chromosome,
                                                      method = "permutation", num.samples = 100,
                                                      use.BLUP = T, model.type = "null")


permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1,
                                      chr = gene_chromosome)

#use.lod = T if want LOD threshold, otherwise I think the number we get is p val?
permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, percentile = 0.85, use.lod = T) # 10.1186/s40168-023-01552-8 uses 85 percentile


# Ensure the directory exists
output_dir <- file.path("Data/Scan_thresholds/Permute_1Gene1Chr", args[3], args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the output file path based on the phenotype of interest
output_threshold <- file.path(output_dir, paste0(as.character(args[1]), "_threshold.rds"))
output_scan <- file.path(output_dir, paste0(as.character(args[1]), "_permuted_scan.rds"))

# Save the threshold and scan as RDS files
saveRDS(permute_threshold, output_threshold)
saveRDS(permuted_scans, output_scan)

