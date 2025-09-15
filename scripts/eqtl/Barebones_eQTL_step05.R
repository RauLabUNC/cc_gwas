#Summary: Each gene from the trait loci has haplotype blocks of eQTL. In each block there are SNPs. Is there at least a significant SNP for this gene that is under the same trait locus?

library(dplyr); library(tidyr)
library(openxlsx)
library(readr)

# ---------------------------------------------------------
# Inputs : genes_mouse.csv, the snp file, the gene names
# Output : a file that has gene names, its eQTL SNPs, is that SNP in the trait locus that gene's in, etc
# ---------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

## 0. paths ------------------------------------------------
#Do this separately for control and ISO groups

#The file with each trait significant loci and genes in there
genes_file <- "Trait_Loci_Genes.xlsx"

#the consolidate file with all snps eqtl for that gene
eqtl_file  <- "control_snps.csv" #This can also be ISO_snps.csv

#outfile name, careful with name
out_file   <- "EQTL_FINAL_RESULTS.csv"

## 1. read tables -----------------------------------------
#The one for eQTL
eqtl = read.csv(eqtl_file, row.names = 1) %>%
  dplyr::select(SNP, pos, chr, trait) %>%
  mutate(gene_sym = trait)

eqtl = data.table(eqtl) #just so later steps will work

#The one for genes
genes <- read.xlsx(genes_file, sheet = "Loci_genes")
genes$gene_sym = genes$Gene

genes = data.table(genes)


## 2. Input the trait locus window of each gene --------------------------
genes[ , `:=`(win_start = locus.Start,
              win_end   = locus.End)]

genes_win <- genes[ , .(gene_sym, Chr, win_start, win_end)]


## 3. join eQTL ↔ gene windows by gene symbol -------------
setkey(eqtl,  gene_sym)
setkey(genes_win, gene_sym)

eqtl_gene <- genes_win[eqtl] #%>%    # keeps all rows of eqtl, joins gene coords
# na.omit()


## 4. compute cis/trans -----------------------------------
eqtl_gene[ , cis_flag :=
             fifelse(
               !is.na(chr) &                        # gene found
                 Chr == chr &                        # same chromosome
                 pos   >= win_start &       # interval overlaps window
                 pos <= win_end,
               "eQTL_here", "no")
] 

final_df = eqtl_gene %>%
  select(-trait)
colnames(final_df) = c("Gene", "gene.Chr", "loci.Start","loci.End", "SNP", "SNP.position", "snp.Chr", "is.this.eQTL.in.loci")


## 5. write result ----------------------------------------
write.csv(final_df, out_file)
#Can also join the two treatment groups files together as a master table for ease of use when assembling evidence for gene selection