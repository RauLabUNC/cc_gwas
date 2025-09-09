# =============================================================================
# Query MouseMine for Gene Phenotype Annotations
#
# Description:
# This script queries the MouseMine database for phenotype annotations of genes
# identified in QTL loci. Uses InterMineR to fetch mammalian phenotype ontology
# annotations for mouse genes.
# =============================================================================

# --- Load Libraries ---
library(InterMineR)
library(data.table)

# --- Load Input Data ---
genes <- read.csv("results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv")
cat("Genes loaded:", nrow(genes), "\n")

# --- Set up MouseMine Connection ---
im <- initInterMine(mine = listMines()["MouseMine"])

# Get available templates
template <- getTemplates(im)

# --- Configure Query for Gene Phenotypes ---
queryGenePath <- getTemplateQuery(
  im = im, 
  name = "_Feature_Phenotype"
)

# Add gene list as query constraint
queryGenePath$where[[4]] <- list(
  path = c("OntologyAnnotation.subject"),
  op = c("LOOKUP"), 
  value = c(paste(genes$mouse_gene_symbol, collapse = ",")), 
  code = c("B")
)

# --- Execute Query ---
cat("Querying MouseMine for phenotypes...\n")
results <- runQuery(im, queryGenePath)
cat("Phenotype annotations retrieved:", nrow(results), "\n")

# --- Format and Save Output ---
# Rename column for consistency
colnames(results)[which(colnames(results) == "OntologyAnnotation.subject.symbol")] <- "mouse_gene_symbol"

# Ensure output directory exists
output_dir <- "data/processed/joinLoci/relational_tables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save phenotype annotations
write.csv(results, file.path(output_dir, "mouseGenePhenotypes.csv"), row.names = FALSE)
cat("Results saved to:", file.path(output_dir, "mouseGenePhenotypes.csv"), "\n")
