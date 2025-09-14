# --- Load Libraries ---
library(InterMineR)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("--input_gene_info"), type = "character", help = "Path to input table of genes in trait loci"),
  make_option(c("--output_ms_gene_phenos"), type = "character", help = "Path to output csv of mouse gene phenotypes")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)

print(opt)
# --- Load Input Data ---
genes <- read.csv(opt$options$input_gene_info)
cat("Genes loaded:", nrow(genes), "\n")
gene_syms <- unique(stats::na.omit(genes$mouse_gene_symbol))
# --- Set up MouseMine Connection ---
im <- initInterMine(mine = listMines()["MouseMine"])

# Get available templates
template <- getTemplates(im)

# --- Configure Query for Gene Phenotypes ---
queryGenePath <- getTemplateQuery(
  im = im, 
  name = "_Feature_Phenotype"
)

# We'll batch the gene list into chunks to avoid long URLs (HTTP 414)
mk_chunk_indices <- function(n, k) ceiling(seq_len(n) / k)
chunk_size <- 200  # adjust if needed
chunks <- split(gene_syms, mk_chunk_indices(length(gene_syms), chunk_size))

# --- Execute Query ---
cat("Querying MouseMine for phenotypes in", length(chunks), "chunks...\n")
results_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  syms <- chunks[[i]]
  qp <- queryGenePath
  qp$where[[4]] <- list(
    path = c("OntologyAnnotation.subject"),
    op = c("LOOKUP"),
    value = c(paste(syms, collapse = ",")),
    code = c("B")
  )
  cat(sprintf("  • Chunk %d/%d (%d genes) ... ", i, length(chunks), length(syms)))
  res <- tryCatch({
    runQuery(im, qp)
  }, error = function(e) {
    warning(sprintf("chunk %d failed: %s", i, e$message))
    data.frame()
  })
  if (NROW(res) > 0) {
    results_list[[i]] <- res
    cat(sprintf("%d rows\n", nrow(res)))
  } else {
    cat("0 rows\n")
  }
  Sys.sleep(0.2)  # be kind to the API
}
results <- data.table::rbindlist(results_list, fill = TRUE)
cat("Phenotype annotations retrieved:", nrow(results), "\n")

# --- Format and Save Output ---
if (nrow(results) > 0) {
  # Rename column for consistency
  colnames(results)[which(colnames(results) == "OntologyAnnotation.subject.symbol")] <- "mouse_gene_symbol"
}

# Save phenotype annotations
write.csv(results, file.path(opt$options$output_ms_gene_phenos), row.names = FALSE)
cat("Results saved to:", file.path(opt$options$output_ms_gene_phenos), "\n")
