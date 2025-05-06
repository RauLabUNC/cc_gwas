library(data.table)
gi <- fread("data/processed/joinLoci/relational_tables/gene_info.csv",
            select = c("mouse_ensembl_id", "mouse_gene_symbol",
                       "description", "gene_biotype",
                       "chromosome_name", "start_position", "end_position",
                       "human_ensembl_id", "human_gene_symbol"))

## 1. genes_mouse -----------------------------------------------------------
genes_mouse <- unique(
  gi[, .(mouse_ensembl_id,
         mouse_gene_symbol,
         description,
         gene_biotype,
         chr = chromosome_name,
         start_bp = start_position,
         end_bp   = end_position)]
)
setorder(genes_mouse, chr, start_bp)
fwrite(genes_mouse, "data/processed/joinLoci/relational_tables/genes_mouse.csv")

## 2. orthology -------------------------------------------------------------
orthology <- unique(
  gi[, .(mouse_ensembl_id,
         human_ensembl_id,
         human_gene_symbol)]
)
fwrite(orthology, "data/processed/joinLoci/relational_tables/orthology.csv")
