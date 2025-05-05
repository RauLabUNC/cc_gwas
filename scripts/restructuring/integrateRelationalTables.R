# integrateRelationalTables.R
# Script to load relational tables, connect eQTL loci, gene info, and disease associations,
# then produce an integrated table for downstream analysis.

library(data.table)

# Directory of relational tables
dir <- "data/processed/joinLoci/relational_tables"

# 1) Load tables
gene_info    <- fread(file.path(dir, "gene_info.csv"))
associations <- fread(file.path(dir, "associations.csv"))
constraints  <- fread(file.path(dir, "constraints.csv"))
tractability<- fread(file.path(dir, "tractability.csv"))
loci         <- fread(file.path(dir, "loci.csv"))

# 2) Inspect basic dimensions and unique counts
cat("Dimensions and unique counts for each table:\n")
for (nm in c("gene_info","associations","constraints","tractability","loci")) {
  dt <- get(nm)
  cat(sprintf("%s: %d rows, %d cols; unique mouse genes: %d; unique human genes: %d\n",
              nm, nrow(dt), ncol(dt),
              ifelse("mouse_ensembl_id" %in% names(dt), length(unique(dt$mouse_ensembl_id)), NA),
              ifelse("human_ensembl_id" %in% names(dt), length(unique(dt$human_ensembl_id)), NA)))
}

# 3) Subset to eQTL loci
eqtl_loci <- loci[Analysis_Type == "eQTL"]
cat(sprintf("Filtered eQTL loci: %d rows, %d unique loci, %d unique traits\n",
            nrow(eqtl_loci), length(unique(eqtl_loci$Locus_ID)), length(unique(eqtl_loci$Dependent_Variable))))

# 4) Merge gene_info with eQTL loci (gene to locus mapping for eQTL)
gene_eqtl <- merge(
  gene_info,
  eqtl_loci[, .(Locus_ID, Dependent_Variable, Chr, Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Peak_SNP_pos_bp, Peak_Significance_Value)],
  by = "Locus_ID",
  all = FALSE
)
cat(sprintf("After merging gene_info with eQTL loci: %d rows\n", nrow(gene_eqtl)))

# 5) Merge with disease associations
gene_eqtl_assoc <- merge(
  gene_eqtl,
  associations,
  by = "human_ensembl_id",
  allow.cartesian = TRUE
)
cat(sprintf("After merging with associations: %d rows; unique diseases: %d\n",
            nrow(gene_eqtl_assoc), length(unique(gene_eqtl_assoc$disease_name))))

# 6) Merge constraints and tractability (allowing cartesian joins)
# Note: constraints and tractability have multiple records per gene, so we allow cartesian
integrated <- merge(
  gene_eqtl_assoc,
  constraints,
  by = "human_ensembl_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)
integrated <- merge(
  integrated,
  tractability,
  by = "human_ensembl_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)
cat(sprintf("Final integrated table: %d rows; columns: %s\n",
            nrow(integrated), paste(names(integrated), collapse = ", "))) 

# 7) Save integrated data
out_file <- file.path(dir, "integrated_eqtl_associations.csv")
fwrite(integrated, out_file)
cat("Saved integrated table to", out_file, "\n")
