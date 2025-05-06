# -----------------------------------------------------------------------------
# Identify Multi‑Trait Candidate Genes with cis‑eQTL Support & NRVM Expression
# -----------------------------------------------------------------------------
#  Selection criteria
#  ------------------
#  1. Mouse gene is protein‑coding.
#  2. Has a mapped human orthologue (non‑empty symbol).
#  3. Gene appears in cis eQTLs.
#  4. Overlaps ≥ 10 miQTL trait loci.
#  5. Occurs in > 20 Trait×Drug combinations **and** avg CPM > 30 in NRVM Ctrl.
# -----------------------------------------------------------------------------
#  Output  : results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv
# -----------------------------------------------------------------------------

# === 0.  Libraries ============================================================
library(data.table)
# library(dplyr)   # Uncomment if you really want tidyverbs

# === 1.  File paths ===========================================================
base   <- "data/processed/joinLoci/relational_tables"
paths  <- list(
  genes_mouse   = fs::path(base, "genes_mouse.csv"),
  orthology     = fs::path(base, "orthology.csv"),
  associations  = fs::path(base, "associations.csv"),
  trait_loci    = fs::path(base, "traitLoci.csv"),
  exp_loci_cis  = fs::path(base, "expLoci_cisflag.csv"),
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
exp_cis       <- read_dt(paths$exp_loci_cis)   # already flagged cis/trans
associations  <- read_dt(paths$associations)  
# === 4.  Filter to protein‑coding genes with human orthologues ===============
genes_pc      <- genes_mouse[gene_biotype == "protein_coding"]
genes_pc_h    <- merge(genes_pc, ortho_mouse2h, by = "mouse_ensembl_id")[human_gene_symbol != ""]

# Keep only genes present in cis eQTLs (dependent variable list) ---------------
genes_pc_h_cis <- genes_pc_h[mouse_gene_symbol %in% exp_cis$Dependent_Variable]

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
GeneSumm <- merge(GeneTrait, NRVM_expr, by = "mouse_gene_symbol", all.x = TRUE)

Final <- GeneSumm[avgenes_cpm > 30][
  order(-n_trait_drug, -avgenes_cpm)]



# === 9. get the pathologies of the top genes 
top_genes <- as.data.frame(Final)[1:10, "mouse_ensembl_id"]

top_genes_human <- ortho_mouse2h |> filter(mouse_ensembl_id %in% top_genes) |> pull(human_gene_symbol)

top_patho <- associations |> filter(symbol %in% top_genes_human)

# === 9.  Write output =========================================================
out_dir <- "results/joinLoci/geneTables"
fs::dir_create(out_dir)
write_dt(Final,
         fs::path(out_dir, "multTrait_cis-eQTL_nrvmExp.csv"))
