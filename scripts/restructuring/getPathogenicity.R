# getGeneInfo.R
# Script to build relational tables of gene information, disease associations,
# genetic constraints, and tractability using Open Targets batch GraphQL queries

# Load required libraries
library(data.table)
library(biomaRt)
library(httr)
library(jsonlite)
library(purrr)
library(dplyr)
library(tibble)

# 1) Read genes-in-loci output
loci_genes <- readRDS("data/processed/joinLoci/relational_tables/genesInLoci.rds")

# 2) Unnest Ensembl IDs per locus
gene_dt <- loci_genes[, .(mouse_ensembl_id = unlist(Ensembl_IDs)), by = .(Locus_ID)]
unique_genes <- unique(gene_dt$mouse_ensembl_id)

# 3) Set up Ensembl mart connection
mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# 4) Query mouse gene info in one go
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

# 5) Query human ortholog mapping in one go
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

# 6) Merge attributes and ortholog mapping, then dedupe
gene_info <- merge(
  mouse_attrs,
  homology,
  by = "mouse_ensembl_id",
  all.x = TRUE
)
gene_info <- unique(gene_info, by = "mouse_ensembl_id")

# 7) Join with locus information
gene_loci_info <- merge(
  gene_dt,
  gene_info,
  by = "mouse_ensembl_id",
  all.x = TRUE
)

# 8) Build GraphQL query for Open Targets
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

# 9) Batch GraphQL requests and parse results
endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
human_ids <- na.omit(unique(gene_info$human_ensembl_id))
# reuse make_batches for GraphQL batches only
graph_batches <- split(human_ids, ceiling(seq_along(human_ids)/100))

assoc_list      <- vector("list", length(graph_batches))
constraint_list <- vector("list", length(graph_batches))
tract_list      <- vector("list", length(graph_batches))


#tract_safe <- tract_list
#tract_list[[i]] <- NA
for (i in seq_along(graph_batches)) {
  if(is.null(tract_list[[i]])){
  ids <- graph_batches[[i]]
  qry <- build_query(ids)
  res <- POST(endpoint, body = list(query = qry), encode = "json")
  stop_for_status(res)
  raw_data <- fromJSON(content(res, as = "text"), flatten = FALSE)$data
  print(i/length(graph_batches))
  # Disease associations
  df_assocs <- imap_dfr(raw_data, ~ {
    rows <- .x$associatedDiseases$rows
    tibble(
      human_ensembl_id  = sub("^g_", "", .y),
      symbol            = .x$approvedSymbol,
      disease_id        = rows$disease$id,
      disease_name      = rows$disease$name,
      association_score = rows$score
    )
  })
  assoc_list[[i]] <- df_assocs
  
  # Genetic constraints
  df_constraints <- map_dfr(raw_data, ~ .x$geneticConstraint, .id = "alias") %>%
    mutate(human_ensembl_id = sub("^g_", "", alias)) %>%
    dplyr::select(-alias)
  constraint_list[[i]] <- df_constraints
  
  # Tractability
  df_tract <- map_dfr(raw_data, ~ .x$tractability, .id = "alias") %>%
    mutate(human_ensembl_id = sub("^g_", "", alias)) %>%
    dplyr::select(-alias)
  tract_list[[i]] <- df_tract
  
  Sys.sleep(1) # polite pause
  if(i %% 20){
    saveRDS(tract_list, "data/processed/joinLoci/relational_tables/geneBatchesAPI/trackList.rds")
  }
  }
}
#tract_list[[153]] <- data.frame()
all_assocs       <- bind_rows(assoc_list)
all_constraints  <- bind_rows(constraint_list)
all_tractability <- bind_rows(tract_list)

# 10) Save outputs
fwrite(gene_loci_info,     "data/processed/joinLoci/relational_tables/gene_info.csv")
fwrite(all_assocs,         "data/processed/joinLoci/relational_tables/associations.csv")
fwrite(all_constraints,    "data/processed/joinLoci/relational_tables/constraints.csv")
fwrite(all_tractability,   "data/processed/joinLoci/relational_tables/tractability.csv")
