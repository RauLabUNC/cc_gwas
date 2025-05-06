library(data.table)

## -------------------------------------------------------------------------
## 1.  LOAD TABLES already on disk
## -------------------------------------------------------------------------
trait   <- fread("data/processed/joinLoci/relational_tables/traitLoci.csv")[Method == "miQTL"]
pos     <- fread("data/processed/joinLoci/relational_tables/pos.csv")
genes   <- fread("data/processed/joinLoci/relational_tables/genes_mouse.csv")[gene_biotype == "protein_coding"]
ortho   <- fread("data/processed/joinLoci/relational_tables/orthology.csv")      # mouse ↔ human

## -------------------------------------------------------------------------
## 2.  Give every locus explicit coordinates
## -------------------------------------------------------------------------
trait_pos <- merge(
  trait[, .(Locus_ID, Position_ID)],   # keep only the bits we need
  pos,                                 # chr, start_bp, end_bp
  by.x = "Position_ID", by.y = "pos_id",
  all.x = TRUE, allow.cartesian = TRUE
)

## -------------------------------------------------------------------------
## 3.  Ensure both tables share the same 3‑column key
## -------------------------------------------------------------------------
## trait_pos already has  chr / start_bp / end_bp
setkey(trait_pos, chr, start_bp, end_bp)

## genes_mouse too
setkey(genes,      chr, start_bp, end_bp)

## -------------------------------------------------------------------------
## 4.  Interval overlap
## -------------------------------------------------------------------------
ov <- foverlaps(
  x      = genes,        # query intervals
  y      = trait_pos,    # subject intervals
  type   = "any",
  nomatch = 0L
)

## ov now has one row per (gene × locus) hit
## -------------------------------------------------------------------------
## 5.  Count distinct loci per gene
## -------------------------------------------------------------------------
gene_hits <- ov[, .(n_trait_loci = uniqueN(Locus_ID)),
                by = .(mouse_ensembl_id, mouse_gene_symbol)]

## 6. keep genes with a human orthologue
gene_hits <- merge(gene_hits, ortho,
                   by = "mouse_ensembl_id", allow.cartesian = TRUE)

setorder(gene_hits, -n_trait_loci)
head(gene_hits, 30)


trait <- fread("data/processed/joinLoci/relational_tables/traitLoci.csv")[Method == "miQTL",
                                                  .(Locus_ID, Dependent_Variable, Treatment)]

## ------------------------------------------------------------------
## 1.  Join ov ↔ trait to bring in Trait & Drug for every hit
## ------------------------------------------------------------------
ov2 <- merge(
  ov,                                # the big gene × locus overlap table
  trait,                             # adds Dependent_Variable + Treatment
  by = "Locus_ID",
  all.x = TRUE, 
  allow.cartesian = T
) |> dplyr::mutate("VarXtreat" = paste0(Dependent_Variable, " - ", Treatment))

## ------------------------------------------------------------------
## 2.  Collapse per gene
## ------------------------------------------------------------------
gene_summ <- ov2[
  , .(
    n_trait_loci = uniqueN(Locus_ID),
    traits       = toString(unique(Dependent_Variable)),   # comma‑sep list
    drugs        = toString(unique(Treatment)),
    VarXtreat    = toString(unique(VarXtreat))
  ),
  by = .(mouse_ensembl_id, mouse_gene_symbol)
]

## ------------------------------------------------------------------
## 3.  Merge back with orthology (if you still need the human IDs)
## ------------------------------------------------------------------
gene_summ <- merge(
  gene_summ,
  orthology,
  by = "mouse_ensembl_id",
  allow.cartesian = TRUE
)

## ------------------------------------------------------------------
## 4.  Inspect the top rows
## ------------------------------------------------------------------
setorder(gene_summ, -n_trait_loci)
gene_summ[1:10]
