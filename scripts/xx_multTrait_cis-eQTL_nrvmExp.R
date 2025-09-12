# -----------------------------------------------------------------------------
# Identify Multi‑Trait Candidate Genes with cis‑eQTL Support & NRVM Expression
# -----------------------------------------------------------------------------
#  Selection criteria
#  ------------------
#  1. Mouse gene is protein‑coding.
#  2. Has a mapped human orthologue (non‑empty symbol).
#  3. avg CPM > 30 in NRVM Ctrl.
# -----------------------------------------------------------------------------
#  Output  : results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv
# -----------------------------------------------------------------------------

# === 0.  Libraries ============================================================
library(data.table)
library(dplyr)   

# === 1.  File paths ===========================================================
base   <- "data/processed/joinLoci/relational_tables"
paths  <- list(
  genes_mouse   = fs::path(base, "genes_mouse.csv"),
  orthology     = fs::path(base, "orthology.csv"),
  associations  = fs::path(base, "associations.csv"),
  trait_loci    = fs::path(base, "traitLoci.csv"),
  # Removed PyLMM dependency - will use miQTL eQTL data instead
  # cis_eQTL_genes= "data/processed/joinLoci/eqtl/PyLMM/VSTdata_Ctrl_and_Iso_2M_CisFlags_250508.csv",
  mouse_pheno   = fs::path(base, "mouseGenePhenotypes.csv"),
  nrvm_counts   = "data/processed/joinLoci/nrvms/bulk_gene.csv",
  nrvm_meta     = "data/processed/joinLoci/nrvms/phenotypes.csv"
)

# === 2.  I/O helpers ==========================================================
read_dt   <- function(x, ...) fread(x, ...)
write_dt  <- function(dt, file) fwrite(dt, file)
message_b <- function(...) cat("\n▸", ..., "\n")

# === 3.  Load core tables =====================================================
genes_mouse   <- read_dt(paths$genes_mouse)
ortho_mouse2h <- read_dt(paths$orthology)
trait_loci    <- read_dt(paths$trait_loci)
associations  <- read_dt(paths$associations)
# Mouse phenotypes may not exist yet - make it optional
if (file.exists(paths$mouse_pheno)) {
  mouse_pheno   <- read_dt(paths$mouse_pheno)
} else {
  message_b("Mouse phenotypes file not found, skipping phenotype annotations")
  mouse_pheno <- data.table()
}  
# Removed PyLMM cis-eQTL dependency - simplified pipeline without cis/trans filtering
# exp_cis       <- read.csv(paths$cis_eQTL_genes, row.names = 1) |> 
#   dplyr::select(
#     mouse_gene_symbol = Gene,
#     Ctrl_cis_flag,
#     Iso_cis_flag
#   )  # already flagged cis/trans
# === 4.  Filter to protein‑coding genes with human orthologues ===============
genes_pc      <- genes_mouse[gene_biotype == "protein_coding"]
genes_pc_h    <- merge(genes_pc, ortho_mouse2h, by = "mouse_ensembl_id")[human_gene_symbol != ""]

# === 5.  Interval overlap: genes vs. miQTL trait loci =========================
# Prepare intervals
trait_int <- trait_loci[Method == "miQTL",
                    .(Locus_ID,
                      Position_ID,
                      trait = Dependent_Variable,
                      drug  = Treatment,
                      chr   = Chr,
                      start = Locus_Start_bp,
                      end   = Locus_End_bp)]
setkey(trait_int, chr, start, end)

gene_int  <- genes_pc_h[, .(mouse_ensembl_id, mouse_gene_symbol,
                            chr, start = start_bp, end = end_bp)]
setkey(gene_int, chr, start, end)

ov <- foverlaps(gene_int, trait_int, type = "any", nomatch = 0L)

# === 6.  Summarise per gene ===================================================
GeneTrait <- ov[, .(
  n_trait_drug = uniqueN(paste(trait, drug, sep = ":")),
  trait_drug   = toString(sort(unique(paste(trait, drug, sep = ":"))))
), by = .(mouse_ensembl_id, mouse_gene_symbol)]


# Function to replace commas with em dashes in a data frame
replace_commas_with_emdashes <- function(df) {
  # Loop through each column in the data frame
  for (col in names(df)) {
    # Check if the column is of character type
    if (is.character(df[[col]])) {
      # Replace commas with em dashes
      df[[col]] <- gsub(",", "—", df[[col]])
    }
  }
  return(df)
}

# Process mouse phenotypes if available
if (nrow(mouse_pheno) > 0) {
  # Apply the function to your data frame
  mouse_pheno_modified <- replace_commas_with_emdashes(mouse_pheno)
  
  mouse_pheno_listed <- mouse_pheno_modified[, .(
    comments = toString(OntologyAnnotation.evidence.comments.description),
    ontology   = toString(OntologyAnnotation.ontologyTerm.name),
    pubmedID = toString(OntologyAnnotation.evidence.publications.pubMedId)
  ), by = .(mouse_gene_symbol)]
} else {
  # Create empty data.table with expected columns
  mouse_pheno_listed <- data.table(
    mouse_gene_symbol = character(),
    comments = character(),
    ontology = character(),
    pubmedID = character()
  )
}

# === 7.  NRVM control expression (CPM) =======================================
NRVM_counts <- read_dt(paths$nrvm_counts)
setnames(NRVM_counts, 1, "mouse_gene_symbol")   # first column is gene symbol
NRVM_meta   <- read_dt(paths$nrvm_meta)

ctrl_ids <- NRVM_meta[treatment == "Ctl", sprintf("BCJ_%s", sample.id)]

# Subset + CPM transform -------------------------------------------------------
counts_mat <- as.matrix(NRVM_counts[, ..ctrl_ids])
lib_sizes  <- colSums(counts_mat)
CPM        <- t( t(counts_mat) / (lib_sizes / 1e6) )  # sweep not needed
avgenes_cpm    <- rowMeans(CPM)

NRVM_expr  <- data.table(mouse_gene_symbol = NRVM_counts$mouse_gene_symbol,
                         avgenes_cpm = avgenes_cpm)

# === 8.  Merge summaries & apply final thresholds ============================
geneSummary <- merge(GeneTrait, NRVM_expr, by = "mouse_gene_symbol", all.x = TRUE) 

geneSummary <- geneSummary[avgenes_cpm > 5][
  order(-n_trait_drug, -avgenes_cpm)]

# === 9. get the pathologies of the top genes 

top_genes_human <- ortho_mouse2h |> filter(mouse_ensembl_id %in% geneSummary$mouse_ensembl_id) |> pull(human_gene_symbol)

top_patho <- associations |> filter(symbol %in% top_genes_human) |> 
  group_by(symbol, human_ensembl_id) |> 
  arrange(desc(association_score)) |> 
  filter(association_score > 0.05) |> 
  as.data.table()

top_patho <- top_patho[, .(
  disease = toString(disease_name)
), by = .(symbol, human_ensembl_id)] |> 
  left_join(ortho_mouse2h)



merged <- geneSummary |> 
  left_join(top_patho, by = "mouse_ensembl_id") |> 
  dplyr::select(-symbol, -human_gene_symbol, -human_ensembl_id) |> 
  dplyr::select(mouse_gene_symbol, n_trait_drug, avgenes_cpm,
                trait_drug, disease, mouse_ensembl_id) |> 
  left_join(mouse_pheno_listed)
  # Removed cis-eQTL join since we're not using PyLMM data
  # left_join(exp_cis)

# === 9.  Write output =========================================================
out_dir <- "results/joinLoci/geneTables"
fs::dir_create(out_dir)
write_dt(merged,
         fs::path(out_dir, "multTrait_cis-eQTL_nrvmExp.csv"))


### Make selected examples 

examples <- merged |> 
  slice_head(n=100) |> 
  group_by(n_trait_drug, trait_drug) |> 
  arrange(desc(avgenes_cpm)) |> 
  slice_head(n = 3) |> 
  ungroup() |> 
  arrange(desc(n_trait_drug), desc(avgenes_cpm))

write.csv(examples, file.path(out_dir, "multTrait_nrvmExp_disease.csv"), row.names = F)

