# Build multi-trait gene summary table with NRVM expression and disease / phenotype annotations
# Refactored from xx_multTrait_cis-eQTL_nrvmExp.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(optparse)
  library(fs)
})

option_list <- list(
  make_option(c("--input_genes_mouse"), type = "character", help = "Path to genes_mouse.csv produced by split species step"),
  make_option(c("--input_orthology"), type = "character", help = "Path to orthology.csv produced by split species step"),
  make_option(c("--input_trait_loci"), type = "character", help = "Path to traitLoci.csv summarizing significant QTL"),
  make_option(c("--input_associations"), type = "character", help = "Path to associations.csv for human gene-disease links"),
  make_option(c("--input_mouse_phenotypes"), type = "character", default = NA, help = "Optional path to mouseGenePhenotypes.csv"),
  make_option(c("--output_gene_table"), type = "character", help = "Output CSV path for full multi-trait gene summary"),
  make_option(c("--output_examples"), type = "character", help = "Output CSV path for top example subset"),
  make_option(c("--cpm_min"), type = "double", default = 5, help = "Minimum average CPM filter (default: 5)"),
  make_option(c("--top_n_examples"), type = "integer", default = 100, help = "Number of top rows to consider before per-stratum selection (default: 100)"),
  make_option(c("--example_per_group"), type = "integer", default = 3, help = "Examples per (n_trait_drug, trait_drug) group (default: 3)")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
print(opt)

message("[multi-trait] Loading inputs ...")

# --- Load core tables ---
read_dt <- function(x, ...) fread(x, ...)
stopifnot(!is.null(opt$input_genes_mouse), !is.null(opt$input_trait_loci))

genes_mouse   <- read_dt(opt$input_genes_mouse)
ortho_mouse2h <- read_dt(opt$input_orthology)
trait_loci    <- read_dt(opt$input_trait_loci)
associations  <- read_dt(opt$input_associations)
mouse_pheno   <- if (!is.na(opt$input_mouse_phenotypes) && file.exists(opt$input_mouse_phenotypes)) read_dt(opt$input_mouse_phenotypes) else data.table()

# --- Protein coding w/ human orthologue ---

genes_pc      <- genes_mouse[gene_biotype == "protein_coding"]
genes_pc_h    <- merge(genes_pc, ortho_mouse2h, by = "mouse_ensembl_id")[human_gene_symbol != ""]

# --- Interval overlap ---
trait_int <- trait_loci[Method == "miQTL",
                        .(Locus_ID, Position_ID, trait = Dependent_Variable, drug = Treatment,
                          chr = Chr, start = Locus_Start_bp, end = Locus_End_bp)]
setkey(trait_int, chr, start, end)

gene_int  <- genes_pc_h[, .(mouse_ensembl_id, mouse_gene_symbol, chr, start = start_bp, end = end_bp)]
setkey(gene_int, chr, start, end)

overlaps <- foverlaps(gene_int, trait_int, type = "any", nomatch = 0L)

GeneTrait <- overlaps[, .(n_trait_drug = uniqueN(paste(trait, drug, sep=":")),
                          trait_drug = toString(sort(unique(paste(trait, drug, sep=":"))))),
                      by = .(mouse_ensembl_id, mouse_gene_symbol)]

# --- Mouse phenotypes summarised (optional) ---
replace_commas_with_emdashes <- function(dt) {
  for (col in names(dt)) if (is.character(dt[[col]])) dt[[col]] <- gsub(",", "—", dt[[col]])
  dt
}

if (nrow(mouse_pheno) > 0) {
  mouse_pheno_mod <- replace_commas_with_emdashes(mouse_pheno)
  mouse_pheno_listed <- mouse_pheno_mod[, .(comments = toString(OntologyAnnotation.evidence.comments.description),
                                            ontology = toString(OntologyAnnotation.ontologyTerm.name),
                                            pubmedID = toString(OntologyAnnotation.evidence.publications.pubMedId)),
                                        by = .(mouse_gene_symbol)]
} else {
  mouse_pheno_listed <- data.table(mouse_gene_symbol=character(), comments=character(), ontology=character(), pubmedID=character())
}

# --- Hardcoded NRVM expression inputs for now ---
NRVM_counts_path <- "data/processed/joinLoci/nrvms/bulk_gene.csv"
NRVM_meta_path   <- "data/processed/joinLoci/nrvms/phenotypes.csv"
NRVM_counts <- read_dt(NRVM_counts_path)
setnames(NRVM_counts, 1, "mouse_gene_symbol")
NRVM_meta   <- read_dt(NRVM_meta_path)
ctrl_ids <- NRVM_meta[treatment == "Ctl", sprintf("BCJ_%s", sample.id)]
counts_mat <- as.matrix(NRVM_counts[, ..ctrl_ids])
lib_sizes  <- colSums(counts_mat)
CPM        <- t( t(counts_mat) / (lib_sizes / 1e6) )
avgenes_cpm <- rowMeans(CPM)
NRVM_expr  <- data.table(mouse_gene_symbol = NRVM_counts$mouse_gene_symbol, avgenes_cpm = avgenes_cpm)

# --- Merge & filter ---
geneSummary <- merge(GeneTrait, NRVM_expr, by = "mouse_gene_symbol", all.x = TRUE)

cpm_min <- opt$cpm_min
geneSummary <- geneSummary[avgenes_cpm > cpm_min][order(-n_trait_drug, -avgenes_cpm)]

# --- Disease associations ---
top_patho <- associations[symbol %in% ortho_mouse2h[mouse_ensembl_id %in% geneSummary$mouse_ensembl_id, human_gene_symbol]]
if (nrow(top_patho) > 0) {
  top_patho <- top_patho[association_score > 0.05]
  if (nrow(top_patho) > 0) {
    top_patho <- top_patho[order(-association_score)]
    top_patho <- top_patho[, .(disease = toString(disease_name)), by = .(symbol, human_ensembl_id)]
    top_patho <- dplyr::left_join(top_patho, ortho_mouse2h,
                                  by = c("symbol" = "human_gene_symbol", "human_ensembl_id" = "human_ensembl_id"))
  } else {
    top_patho <- data.table(mouse_ensembl_id = character(), disease = character())
  }
} else {
  top_patho <- data.table(mouse_ensembl_id = character(), disease = character())
}

merged <- dplyr::left_join(geneSummary, top_patho, by = "mouse_ensembl_id") |> 
  dplyr::select(-dplyr::any_of(c("symbol", "human_gene_symbol", "human_ensembl_id"))) |> 
  dplyr::select(mouse_gene_symbol, n_trait_drug, avgenes_cpm, trait_drug, disease, mouse_ensembl_id) |> 
  dplyr::left_join(mouse_pheno_listed, by = "mouse_gene_symbol")

# --- Write outputs ---
out_full <- opt$output_gene_table
out_examples <- opt$output_examples
if (is.null(out_full) || is.null(out_examples)) stop("Both --output_gene_table and --output_examples must be provided.")
dir_create(path_dir(out_full))

fwrite(merged, out_full)

# Examples subset
examples <- merged |> 
  dplyr::slice_head(n = opt$top_n_examples) |> 
  dplyr::group_by(n_trait_drug, trait_drug) |> 
  dplyr::arrange(dplyr::desc(avgenes_cpm)) |> 
  dplyr::slice_head(n = opt$example_per_group) |> 
  dplyr::ungroup() |> 
  dplyr::arrange(dplyr::desc(n_trait_drug), dplyr::desc(avgenes_cpm))

fwrite(examples, out_examples)

message(sprintf("[multi-trait] Wrote full table: %s (%d rows)", out_full, nrow(merged)))
message(sprintf("[multi-trait] Wrote examples: %s (%d rows)", out_examples, nrow(examples)))
