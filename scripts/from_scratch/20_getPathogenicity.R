# =============================================================================
# Fetch Gene Information and Disease Associations from Open Targets
#
# Description:
# This script builds relational tables of gene information, disease associations,
# genetic constraints, and drug tractability. It queries BioMart for ortholog
# mapping and Open Targets GraphQL API for comprehensive gene annotations.
# =============================================================================

# --- Load Libraries ---
library(data.table)
library(biomaRt)
library(httr)
library(jsonlite)
library(purrr)
library(dplyr)
library(tibble)

# --- Read Genes in Loci ---
loci_genes <- readRDS("data/processed/joinLoci/relational_tables/genesInLoci.rds")

# Unnest Ensembl IDs per locus
gene_dt <- loci_genes[, .(mouse_ensembl_id = unlist(Ensembl_IDs)), by = .(pos_id)]
unique_genes <- unique(gene_dt$mouse_ensembl_id)
cat("Unique genes to process:", length(unique_genes), "\n")

# --- Query Mouse Gene Information from BioMart ---
mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Query mouse gene attributes
mouse_attrs <- getBM(
  attributes = c(
    "ensembl_gene_id", "external_gene_name", "description",
    "gene_biotype", "chromosome_name",
    "start_position", "end_position"
  ),
  filters = "ensembl_gene_id",
  values = unique_genes,
  mart = mouse_mart
)
setDT(mouse_attrs)
setnames(mouse_attrs,
         old = c("ensembl_gene_id","external_gene_name"),
         new = c("mouse_ensembl_id","mouse_gene_symbol"))

# --- Query Human Ortholog Mapping ---
homology <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name"
  ),
  filters = "ensembl_gene_id",
  values = unique_genes,
  mart = mouse_mart
)
setDT(homology)
setnames(homology,
         old = c("ensembl_gene_id",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name"),
         new = c("mouse_ensembl_id","human_ensembl_id","human_gene_symbol"))

# --- Merge and Deduplicate Gene Information ---
gene_info <- merge(
  mouse_attrs,
  homology,
  by = "mouse_ensembl_id",
  all.x = TRUE
)
gene_info <- unique(gene_info, by = "mouse_ensembl_id")

# --- Join with Locus Information ---
gene_loci_info <- merge(
  gene_dt,
  gene_info,
  by = "mouse_ensembl_id",
  all.x = TRUE
)

# --- Build GraphQL Query Function for Open Targets ---
build_query <- function(ids) {
  blocks <- map_chr(ids, ~ {
    alias <- paste0("g_", gsub("-", "_", .x))
    sprintf(
      "%s: target(ensemblId: \"%s\") {\n  id\n  approvedSymbol\n  biotype\n  geneticConstraint { constraintType exp obs score oe oeLower oeUpper }\n  tractability { label modality value }\n  associatedDiseases(page: { size: 100, index: 0 }) { rows { disease { id name } score } }\n}",
      alias, .x
    )
  })
  paste0("query multiTarget {\n", paste(blocks, collapse = "\n"), "\n}")
}

# --- Batch Process Open Targets API Queries ---
endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
human_ids <- na.omit(unique(gene_info$human_ensembl_id))
cat("Human genes to query:", length(human_ids), "\n")

# Split into batches of 100 for API limits
graph_batches <- split(human_ids, ceiling(seq_along(human_ids) / 100))
cat("Processing in", length(graph_batches), "batches\n")

# Initialize result lists
assoc_list <- vector("list", length(graph_batches))
constraint_list <- vector("list", length(graph_batches))
tract_list <- vector("list", length(graph_batches))

# Process each batch
for (i in seq_along(graph_batches)) {
  cat("Processing batch", i, "of", length(graph_batches), "\n")
  
  ids <- graph_batches[[i]]
  qry <- build_query(ids)
  
  # Make API request with error handling
  tryCatch({
    res <- POST(endpoint, body = list(query = qry), encode = "json")
    stop_for_status(res)
    raw_data <- fromJSON(content(res, as = "text"), flatten = FALSE)$data
    
    # Parse disease associations
    df_assocs <- imap_dfr(raw_data, ~ {
      rows <- .x$associatedDiseases$rows
      tibble(
        human_ensembl_id = sub("^g_", "", .y),
        symbol = .x$approvedSymbol,
        disease_id = rows$disease$id,
        disease_name = rows$disease$name,
        association_score = rows$score
      )
    })
    assoc_list[[i]] <- df_assocs
    
    # Parse genetic constraints
    df_constraints <- map_dfr(raw_data, ~ .x$geneticConstraint, .id = "alias") |>
      mutate(human_ensembl_id = sub("^g_", "", alias)) |>
      dplyr::select(-alias)
    constraint_list[[i]] <- df_constraints
    
    # Parse tractability
    df_tract <- map_dfr(raw_data, ~ .x$tractability, .id = "alias") |>
      mutate(human_ensembl_id = sub("^g_", "", alias)) |>
      dplyr::select(-alias)
    tract_list[[i]] <- df_tract
    
  }, error = function(e) {
    cat("Error in batch", i, ":", e$message, "\n")
    assoc_list[[i]] <- data.frame()
    constraint_list[[i]] <- data.frame()
    tract_list[[i]] <- data.frame()
  })
  
  # Rate limiting
  Sys.sleep(1)
  
  # Periodic checkpoint save (every 20 batches)
  if (i %% 20 == 0) {
    checkpoint_dir <- "data/processed/joinLoci/relational_tables/geneBatchesAPI"
    if (!dir.exists(checkpoint_dir)) {
      dir.create(checkpoint_dir, recursive = TRUE)
    }
    saveRDS(list(assoc = assoc_list, constraint = constraint_list, tract = tract_list),
            file.path(checkpoint_dir, paste0("checkpoint_batch_", i, ".rds")))
    cat("  Checkpoint saved at batch", i, "\n")
  }
}

# Combine all results
cat("Combining results...\n")
all_assocs <- bind_rows(assoc_list)
all_constraints <- bind_rows(constraint_list)
all_tractability <- bind_rows(tract_list)

# --- Save Outputs ---
# Ensure output directory exists
output_dir <- "data/processed/joinLoci/relational_tables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

fwrite(gene_loci_info, file.path(output_dir, "gene_info.csv"))
fwrite(all_assocs, file.path(output_dir, "associations.csv"))
fwrite(all_constraints, file.path(output_dir, "constraints.csv"))
fwrite(all_tractability, file.path(output_dir, "tractability.csv"))

cat("\nFiles saved successfully:\n")
cat("  - gene_info.csv:", nrow(gene_loci_info), "genes\n")
cat("  - associations.csv:", nrow(all_assocs), "associations\n")
cat("  - constraints.csv:", nrow(all_constraints), "constraints\n")
cat("  - tractability.csv:", nrow(all_tractability), "tractability records\n")
