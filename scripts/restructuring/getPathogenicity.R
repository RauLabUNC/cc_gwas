# =============================================================================
# Script to query MyGene.info for human orthologs of a mouse QTL region,
# then fetch gene-level data from Open Targets via GraphQL.
# =============================================================================

# --- Setup -------------------------------------------------------------------
# Install dependencies if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("mygene", quietly = TRUE)) BiocManager::install("mygene")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
install_if_missing(c("dplyr", "httr", "jsonlite", "purrr"))

# Load libraries
library(dplyr)
library(httr)
library(jsonlite)
library(purrr)
library(mygene)

# --- 1) Read and extract region ----------------------------------------------
sig_regions <- read.csv("results/sig_regions/all_significant_regions_summary.csv")

# Choose the 5th entry (modify index as needed)
region <- sig_regions %>% slice(n = 5)

chr   <- as.character(region$chr)
start <- region$upper_pos_lod_drop * 1e6
end   <- region$lower_pos_lod_drop * 1e6
interval <- sprintf("chr%s:%d-%d", chr, start, end)

# --- 2) Query MyGene.info for mouse genes in the region ---------------------
mouse_query <- GET(
  url = "http://mygene.info/v3/query",
  query = list(
    q          = interval,
    species    = "mouse",
    fields     = "symbol,agr.orthologs,ensembl.gene",
    size       = 1000,
    ensemblonly = TRUE
  )
)
stop_for_status(mouse_query)
mouse_hits <- fromJSON(content(mouse_query, as = "text"), flatten = TRUE)$hits

# --- 3) Extract best human ortholog Entrez IDs ------------------------------
orthologs <- mouse_hits %>%
  dplyr::select(mouse_symbol = symbol, agr.orthologs) %>%
  unnest(agr.orthologs) %>%
  filter(taxid == 9606) %>%
  group_by(mouse_symbol) %>%
  slice_max(algorithmsmatch, n = 1) %>%
  ungroup() %>%
  transmute(
    mouse_symbol,
    human_entrez = gene_id
  )

# --- 4) Map human Entrez IDs to Ensembl IDs -------------------------------
human_info <- getGenes(
  geneids = orthologs$human_entrez,
  scopes  = "entrezgene",
  species = "human",
  fields  = "ensembl.gene",
  as      = "list"
)
# Build tibble of Entrez → Ensembl
human_map <- tibble(
  entrez = names(human_info),
  info   = human_info
) %>%
  mutate( human_ensembl = map_chr(info, ~ if (!is.null(.x$ensembl.gene)) .x$ensembl.gene$gene else NA_character_)) %>%
  select(entrez, human_ensembl)

# Join back and drop NAs
orthologs <- orthologs %>%
  left_join(human_map, by = c("human_entrez" = "entrez")) %>%
  filter(!is.na(human_ensembl))
human_ids <- orthologs$human_ensembl

# --- 5) Build GraphQL query for Open Targets -------------------------------
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
query_string <- build_query(human_ids[1:100])  # adjust slice as needed

# --- 6) POST to Open Targets GraphQL API ------------------------------------
endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
resp     <- POST(endpoint, body = list(query = query_string), encode = "json")
stop_for_status(resp)

# --- 7) Parse and tidy results ---------------------------------------------
raw_data <- fromJSON(content(resp, as = "text"), flatten = FALSE)$data

# Disease associations
df_assocs <- imap_dfr(raw_data, ~ {
  df <- .x$associatedDiseases$rows
  tibble(
    ensemblId   = sub("^g_", "", .y),
    symbol      = .x$approvedSymbol,
    diseaseId   = df$disease$id,
    diseaseName = df$disease$name,
    assocScore  = df$score
  )
})

# Genetic constraints
df_constraints <- map_dfr(raw_data, ~ .x$geneticConstraint, .id = "alias") %>%
  mutate(ensemblId = sub("^g_", "", alias)) %>%
  select(-alias)

# Tractability
df_tractability <- map_dfr(raw_data, ~ .x$tractability, .id = "alias") %>%
  mutate(ensemblId = sub("^g_", "", alias)) %>%
  select(-alias)

# --- 8) Return or save outputs ---------------------------------------------
list(
  orthologs     = orthologs,
  associations  = df_assocs,
  constraints   = df_constraints,
  tractability  = df_tractability
)
