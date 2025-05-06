library(data.table)
library(dplyr)
library(biomaRt)

# 1) Read & prep loci table ---------------------------------------------------
pos_dt <- fread(
  "data/processed/joinLoci/relational_tables/pos.csv",
  colClasses = list(
    character = "pos_id",
    character = "chr",
    integer   = c("start_bp","end_bp")
  )
) |> 
  drop_na()
# rename for foverlaps
setnames(pos_dt,
         old = c("chr","start_bp","end_bp"),
         new = c("chromosome_name","start","end"))
# set key for overlap
setkey(pos_dt, chromosome_name, start, end)

# 2) Fetch genes & filter ------------------------------------------------------
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)
genes_df <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  mart = mart
)
genes_dt <- as.data.table(genes_df)[
  !grepl("^Gm", external_gene_name) & !grepl("Rik$", external_gene_name) & external_gene_name != "",
  .(
    chromosome_name,
    start = start_position,
    end   = end_position,
    ensembl_gene_id,
    gene_symbol = external_gene_name
  )
]
setkey(genes_dt, chromosome_name, start, end)

# 3) Overlap loci and genes ----------------------------------------------------
# foverlaps requires the first table to have start/end
overlap_dt <- foverlaps(
  genes_dt,
  pos_dt,
  nomatch = 0L
)

# 4) Build list per locus ------------------------------------------------------
# For each locus, collect a vector of unique gene symbols and Ensembl IDs
locus_list_dt <- overlap_dt[
  , .(
    Gene_Symbols   = list(unique(gene_symbol)),
    Ensembl_IDs    = list(unique(ensembl_gene_id))
  ), by = pos_id
]

saveRDS(locus_list_dt, "data/processed/joinLoci/relational_tables/genesInLoci.rds")
