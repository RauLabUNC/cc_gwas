library(dplyr)
library(writexl)
library(openxlsx)
library(stringr)
library(tibble)

##ANH'S NOT TO BRIAN: I did 3 rounds of edits, the final packet that our lab used was 250605. I condensed it to a block of code as if I did it just on=ce

#To edit folder, make a copy of the original one and name it sth new, then work on that folder only
# Edits = add CPM, log2fc, pval, cis tag, #strains gene is expressed in, protein info about cellular location and known protein function

#Info to add/edit---------------------------------------------------------------
#CPM info
ctrl_cpm = read.csv("Data/Processed/Ctrl_CPM.csv", row.names = 1, check.names = F)
iso_cpm = read.csv("Data/Processed/Iso_CPM.csv", row.names = 1, check.names = F)

gene_ctrl = apply(ctrl_cpm, MARGIN=2, mean) %>% as.data.frame() %>%
  rownames_to_column("Gene")
colnames(gene_ctrl)[2] = "Ctrl.Ave.CPM"

gene_iso = apply(iso_cpm, MARGIN=2, mean) %>% as.data.frame() %>%
  rownames_to_column("Gene")
colnames(gene_iso)[2] = "Iso.Ave.CPM"

#DESeq info
deseq = read.csv("Results/DESeq2/deseq2_13kGenes_SexCovar_250521.csv", row.names = 1) #This file is actually in /proj/raulab/users/Anh/CC-miQTL-clean/Results/DESeq2
colnames(deseq)[3] = "deseq.pval"

#Cis eqtl tag
cis_eQTL = read.csv("/proj/raulab/users/Anh/CC-miQTL_VST/Results/Genes_with_cis_eQTL/in_TraitLocus/VST_1Gene1Chr_both_drugs_eqTL_in_locus.csv", row.names = 1)

#Number of strains expressing the gene
num_strains = read.csv("Data/Processed/genes_expressed_in_2plus_strains.csv", row.names = 1, check.names = F)
colnames(num_strains) = c("Gene","Num.Strains.Expressed")

#Got this table from UniProt, all proteins in mouse with protein existence evidence
protein_info = read.delim("Data/Raw/uniprotkb_existence_1_AND_taxonomy_id_1_2025_05_28.tsv", check.names = F) %>%
  mutate(prot_id = stringr::str_remove(STRING, ";.*"),
         STRING.url = if_else(prot_id=="",
                              NA_character_,
                              paste0("https://string-db.org/network/", prot_id)
         )
  )
colnames(protein_info)[4:6] = c("Protein.function", "Subcellular.location", "Tissue.specificity")

##THEN WE CAN JUST JOIN THE COLUMNS AND ARRANGE IN THE ORDER WE WANT



#Input the right file-----------------------------------------------------------

#Pull the files from Brian's packet
loci_folders = list.files("WHERE BRIAN'S FILES ARE")


for (i in c(1:8, 10:12)){ #skip cluster #9 because there's no sig loci genes we care
  
  locus_folder = loci_folders[i]
  loci_name = str_extract(locus_folder, "(?<=locus_).*")
  
  cluster_info_file = paste0("WHERE BRIAN'S FILES ARE", locus_folder, "/gene_info_cluster_", loci_name, ".xlsx")
  
  
  #The current table
  brian_excel = loadWorkbook(cluster_info_file)
  
  current_info = read.xlsx(cluster_info_file, sheet = "AllGenesInCluster") %>%
    select(-Ave_Exp_Ctrl_F, -Ave_Exp_Ctrl_M, -Ave_Exp_Iso_F, -Ave_Exp_Iso_M)
  
  Update_info = current_info %>%
    left_join(x=., gene_ctrl, join_by(Mouse.Gene.Symbol == Gene)) %>%
    left_join(x=., gene_iso, join_by(Mouse.Gene.Symbol == Gene)) %>%
    left_join(x=., deseq, join_by(Mouse.Gene.Symbol == Gene)) %>%
    left_join(x=., cis_eQTL, join_by(Mouse.Gene.Symbol == Gene)) %>%
    left_join(x=., num_strains, join_by(Mouse.Gene.Symbol == Gene)) %>%
    relocate(c(Ctrl.Ave.CPM, Iso.Ave.CPM, log2Fold, deseq.pval , Ctrl_eQTL_here, Iso_eQTL_here), .after = "Avg.NRVM.CPM.(Ctrl)")

  
  Update_info_newSheet = Update_info %>%
    left_join(x=., y=protein_info, join_by(Mouse.Gene.Symbol == 'Gene Names (primary)')) %>%
    dplyr::select(Mouse.Gene.Symbol, Protein.function, Subcellular.location, Tissue.specificity, STRING.url) %>%
    distinct() %>%
    filter(! Protein.function == "") %>%
    filter(rowSums(!is.na(.[, -1])) > 0)
  
  
  
  # Write back the modified sheet
  writeData(brian_excel, sheet = "AllGenesInCluster", x = Update_info)
  addWorksheet(brian_excel, "ProteinInfo")
  writeData(brian_excel, sheet = "ProteinInfo", x = Update_info_newSheet)
  
  # freeze the top header row so it always remains visible
  freezePane(
    brian_excel,
    sheet           = "AllGenesInCluster",
    firstActiveRow  = 2
  )
  
  
  # wrap all text in every cell, so long text columns wrap instead of forcing huge row heights
  wrapStyle <- createStyle(wrapText = TRUE)
  addStyle(
    brian_excel,
    sheet = "AllGenesInCluster",
    style = wrapStyle,
    rows  = 1:(nrow(Update_info)),
    cols  = 1:ncol(Update_info),
    gridExpand = TRUE
  )
  
  # auto-fit every column width (approximate)
  setColWidths(
    brian_excel,
    sheet  = "AllGenesInCluster",
    cols  = 1:ncol(Update_info),
    widths = "auto"
  )
  
  
  # Save workbook without deleting other sheets
  saveWorkbook(brian_excel, cluster_info_file, overwrite = TRUE)
  
}



