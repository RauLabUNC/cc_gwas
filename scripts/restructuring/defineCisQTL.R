#!/usr/bin/env Rscript
# ---------------------------------------------------------
# Tag eQTL loci as cis / trans (0.5 Mb window)
# Inputs : genes_mouse.csv, expLoci.csv
# Output : expLoci_cisflag.csv (same as expLoci + cis_flag column)
# ---------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

## 0. paths ------------------------------------------------
genes_file <- "data/processed/joinLoci/relational_tables/genes_mouse.csv"
eqtl_file  <- "data/processed/joinLoci/relational_tables/expLoci.csv"
out_file   <- "data/processed/joinLoci/relational_tables/expLoci_cisflag.csv"

## 1. read tables -----------------------------------------
genes <- fread(genes_file,      # 25 k rows
               select = c("mouse_gene_symbol", "chr",
                          "start_bp", "end_bp"))

eqtl  <- fread(eqtl_file,       # 1.4 M rows
               select = c("Locus_ID", "Position_ID", "Dependent_Variable", "Chr",
                          "Locus_Start_bp", "Locus_End_bp"))
# standardise gene symbols to upper‑case
genes[ , gene_sym := toupper(mouse_gene_symbol)]
eqtl[  , gene_sym := toupper(Dependent_Variable)]

## 2. create ±0.5 Mb windows around each gene -------------
flank <- 5e5  # 0.5 Mb
genes[ , `:=`(win_start = pmax(start_bp - flank, 0L),
              win_end   = end_bp   + flank)]

genes_win <- genes[ , .(gene_sym, chr, win_start, win_end)]

## 3. join eQTL ↔ gene windows by gene symbol -------------
setkey(eqtl,   gene_sym)
setkey(genes_win, gene_sym)

eqtl_gene <- genes_win[eqtl]    # keeps all rows of eqtl, joins gene coords

## 4. compute cis/trans -----------------------------------
eqtl_gene[ , cis_flag :=
             fifelse(
               !is.na(chr) &                        # gene found
                 Chr == chr &                        # same chromosome
                 Locus_End_bp   >= win_start &       # interval overlaps window
                 Locus_Start_bp <= win_end,
               "cis", "trans")
] 

eqtl_gene <- eqtl_gene |> filter(cis_flag == "cis")


## 5. write result ----------------------------------------
fwrite(eqtl_gene[
  , .(Locus_ID, Position_ID, Method = NULL, Analysis_Type = NULL,
      Dependent_Variable, Chr,
      Locus_Start_bp, Locus_End_bp, cis_flag)],
  out_file)

cat("Wrote", out_file, "with", nrow(eqtl_gene), "rows\n")
