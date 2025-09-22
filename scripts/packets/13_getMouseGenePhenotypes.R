# --- Load Libraries ---
#library(InterMineR)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("--input_gene_info"), type = "character", help = "Path to input table of genes in trait loci"),
  make_option(c("--output_ms_gene_phenos"), type = "character", help = "Path to output csv of mouse gene phenotypes")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)

print(opt)
##### Dummy section to just load in old phenos and skip querying #####
results <- read.csv("/proj/raulab/projects/cc_gwas/data/processed/joinLoci/relational_tables/mouseGenePhenotypes.csv")

# --- Format and Save Output ---
if (nrow(results) > 0) {
  # Rename column for consistency
  colnames(results)[which(colnames(results) == "OntologyAnnotation.subject.symbol")] <- "mouse_gene_symbol"
}

# Save phenotype annotations
write.csv(results, file.path(opt$options$output_ms_gene_phenos), row.names = FALSE)
cat("Results saved to:", file.path(opt$options$output_ms_gene_phenos), "\n")
