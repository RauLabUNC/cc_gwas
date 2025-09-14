#Summary: Each gene from the trait loci has haplotype blocks of eQTL. In each block there are SNPs. Is there at least a significant SNP for this gene that is under the same trait locus?

library(dplyr); library(tidyr)
library(openxlsx)
library(readr)

# ---------------------------------------------------------
# Tag eQTL loci as cis / trans (1 Mb window)
# Inputs : genes_mouse.csv, the snp file, the gene names
# Output : a file that has gene names, its eQTL SNPs, is that SNP in the trait locus that gene's in, etc
# ---------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

## 0. paths ------------------------------------------------
  #Do the same thing for iso data, make sure you have the right name

#The file with each brian's trait loci and genes in there. I made it manually
genes_file <- "Data/Processed/Trait_Loci_Genes.xlsx"

#the file with all snps eqtl for that gene
eqtl_file  <- "Results/Consolidate_SNPs/GeneOwnChr/250602-VST_614/control_snps.csv"

#outfile name, careful with name
out_file   <- "Results/Genes_with_cis_eQTL/in_TraitLocus/VST_1Gene1Chr_control_eqTL.csv"

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

#eqtl_gene <- eqtl_gene |> filter(cis_flag == "cis")


## 5. write result ----------------------------------------
write.csv(final_df, out_file)



################################################################################
#Use another script to merge ctrl and iso into 1 table, I had to run this script as its own job bc R studio crashes. But here's the code

#BY TRAIT LOCI only did this on 614 genes---------------------------------------
ctrl_file = read.csv(file = "Results/Genes_with_cis_eQTL/in_TraitLocus/VST_1Gene1Chr_control_eqTL.csv", row.names = 1, check.names = F) %>%
  select(Gene, is.this.eQTL.in.loci) %>%
  group_by(Gene) %>%
  summarise(is.this.eQTL.in.loci2 = ifelse(any(is.this.eQTL.in.loci == "eQTL_here"), "eQTL_here", "no"))
colnames(ctrl_file) = c("Gene","Ctrl_eQTL_here")


iso_file = read.csv(file = "Results/Genes_with_cis_eQTL/in_TraitLocus/VST_1Gene1Chr_iso_eqTL.csv", row.names = 1, check.names = F) %>%
  select(Gene, is.this.eQTL.in.loci) %>%
  group_by(Gene) %>%
  summarise(is.this.eQTL.in.loci = ifelse(any(is.this.eQTL.in.loci == "eQTL_here"), "eQTL_here", "no"))
colnames(iso_file) = c("Gene","Iso_eQTL_here")


Combined = full_join(ctrl_file, iso_file, join_by(Gene))%>% 
  mutate(across(everything(), ~replace_na(., "no")))

write.csv(Combined, "Results/Genes_with_cis_eQTL/in_TraitLocus/VST_1Gene1Chr_both_drugs_eqTL_in_locus.csv")


