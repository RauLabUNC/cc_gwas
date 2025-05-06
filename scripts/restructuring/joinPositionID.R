# this script is mean to consolidate all positional IDs from the integrated 
# ... cc gwas study

library(data.table)

## 1. Grab Position_ID from both locus tables -------------------------------
trait_pos <- fread(
  "data/processed/joinLoci/relational_tables/traitLoci.csv",
  select = "Position_ID"               # ← read only the needed column
)
eqtl_pos  <- fread(
  "data/processed/joinLoci/relational_tables/expLoci.csv",
  select = "Position_ID"
)

## 2. Union‑deduplicate & parse chr:start‑end -------------------------------
pos <- unique(rbindlist(list(trait_pos, eqtl_pos)))

# Split "chr:start-end" → chr, start, end
pos[, c("chr","start_bp","end_bp") := 
      tstrsplit(Position_ID, "[:-]", type.convert=TRUE, keep=1:3)]

# Optional extras
pos[, length_bp := end_bp - start_bp + 1L]
setorder(pos, chr, start_bp)

## 3. Rename & write --------------------------------------------------------
setnames(pos, "Position_ID", "pos_id")
fwrite(pos, "data/processed/joinLoci/relational_tables/pos.csv")
