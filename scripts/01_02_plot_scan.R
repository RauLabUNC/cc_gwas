library(miqtl)

# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
args <- c("Fibroblast", "iso")
phenotype_of_interest <- args[1]


# Load the data
data_file <- file.path("data/processed/scans", args[2], paste0(as.character(phenotype_of_interest), "_scan_results.rds"))
data <- readRDS(data_file)

# Plot the results
output_dir <- file.path("results/genome_scans", args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

png_name <- file.path(output_dir, paste0(as.character(phenotype_of_interest), "_scan_results.png"))
png(file = png_name,
    width = 8, 
    height = 4,
    units = "in",
    res = 600)

genome.plotter.whole(scan.list=list(ROP = data))

dev.off()
