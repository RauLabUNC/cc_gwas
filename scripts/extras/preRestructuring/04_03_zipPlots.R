# Load necessary library
library(zip)

# Define the directory containing the plots
plot_dir <- "results/locusZoomPlots"

# List all files in the subdirectories, excluding those that end in dummy.txt
plot_files <- list.files(plot_dir, recursive = TRUE, full.names = TRUE)
plot_files <- plot_files[!grepl("dummy\\.txt$", plot_files)]

# Define the output zip file
zip_file <- file.path(plot_dir, "allPlots.zip")

# Create the zip file containing all the plot files
zip(zipfile = zip_file, files = plot_files)

# Print a message indicating completion
cat("All plot files have been zipped into", zip_file, "\n")