library(data.table)
library(dplyr)
library(biomaRt)

# 1) Read & prep loci table ---------------------------------------------------
loci_dt <- fread(
  "data/processed/joinLoci/relational_tables/loci.csv",
  colClasses = list(
    character = "Locus_ID",
    character = "Chr",
    integer   = c("Locus_Start_bp","Locus_End_bp")
  )
)

loci_deDup <- loci_dt |> group_by(Chr, Locus_Start_bp, Locus_End_bp) |> slice_head(n = 1)  |> as.data.table()
# rename for foverlaps
setnames(loci_deDup,
         old = c("Chr","Locus_Start_bp","Locus_End_bp"),
         new = c("chromosome_name","start","end"))
# set key for overlap
setkey(loci_deDup, chromosome_name, start, end)

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
  loci_deDup,
  nomatch = 0L
)

# 4) Build list per locus ------------------------------------------------------
# For each locus, collect a vector of unique gene symbols and Ensembl IDs
locus_list_dt <- overlap_dt[
  , .(
    Gene_Symbols   = list(unique(gene_symbol)),
    Ensembl_IDs    = list(unique(ensembl_gene_id))
  ), by = Locus_ID
]

locus_list_dt[ , locus_range := tstrsplit(Locus_ID, "_", keep=1) ]

saveRDS(locus_list_dt, "data/processed/joinLoci/relational_tables/genesInLoci.rds")
