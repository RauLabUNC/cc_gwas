# --- 4. Define Table Generation Functions --- ####

# Function to generate Gene Information Excel File with multiple sheets
generate_gene_info_excel <- function(genes_in_locus_dt, 
                                     loci_info_df, # Changed: Now expects a data frame of all loci in the cluster
                                     output_dir,
                                     full_mouse_pheno_data = mouse_pheno,
                                     full_human_assoc_data = associations,
                                     orthology_data = ortho_mouse2h) {
  
  if (nrow(genes_in_locus_dt) == 0) {
    message("No genes found in the specified locus for Excel generation.")
    return(NULL)
  }
  
  # --- Create the main summary table with locus presence columns ---
  ensembl_ids_in_locus <- genes_in_locus_dt$mouse_ensembl_id
  locus_gene_summary <- merged_gene_info[mouse_ensembl_id %in% ensembl_ids_in_locus, ]
  
  # Get founder mutations for genes in locus
  founder_muts_gene_level <- muNoSplice[Gene %in% locus_gene_summary$mouse_gene_symbol, 
                                        .(mouse_gene_symbol = Gene, founder_gene_mutations = Mutations)]
  
  if (nrow(founder_muts_gene_level) > 0) {
    locus_gene_summary <- merge(locus_gene_summary, founder_muts_gene_level, by = "mouse_gene_symbol", all.x = TRUE)
  } else {
    locus_gene_summary[, founder_gene_mutations := NA_character_]
  }
  
  # Create a data table with main gene information
  main_genes_table <- locus_gene_summary[, .(
    `Mouse Gene Symbol` = mouse_gene_symbol,
    `Mouse Ensembl ID` = mouse_ensembl_id,
    `# Traits Associated (miQTL)` = n_trait_drug,
    `Associated Traits (Drug)` = trait_drug,
    `Avg NRVM CPM (Ctrl)` = round(avgenes_cpm, 2),
    `Human Disease Summary` = disease, # This is the summary string from merged_gene_info
    `Mouse Phenotype Summary (MGI)` = ontology, # Summary string
    `Supporting Publications (MGI)` = pubmedID,
    `Founder Gene Mutations (CC)` = founder_gene_mutations
  )]
  
  # Add gene coordinates
  gene_coords <- genes_mouse[mouse_ensembl_id %in% main_genes_table$`Mouse Ensembl ID`, 
                             .(mouse_ensembl_id, chr, start_bp, end_bp)]
  setnames(gene_coords, "mouse_ensembl_id", "Mouse Ensembl ID")
  main_genes_table <- merge(main_genes_table, gene_coords, by = "Mouse Ensembl ID", all.x = TRUE)
  
  # Set initial column order
  initial_cols <- c("Mouse Gene Symbol", "Mouse Ensembl ID", "chr", "start_bp", "end_bp", 
                    "Avg NRVM CPM (Ctrl)", "# Traits Associated (miQTL)")
  
  # --- Add locus presence columns (T/F) for each locus in the cluster ---
  # Group loci by trait and drug for unique identifiers
  loci_info_df$locus_id <- paste0(loci_info_df$trait, "_", loci_info_df$drug, 
                                  "_chr", loci_info_df$chr, "_", round(loci_info_df$peak_pos, 2), "Mb")
  
  # For each locus, create a column indicating if each gene is present in that locus
  for (i in 1:nrow(loci_info_df)) {
    current_locus <- loci_info_df[i, ]
    locus_col_name <- paste0("In_", current_locus$locus_id)
    
    # Get the genes in this specific locus
    locus_start_bp <- min(current_locus$upper_pos_lod_drop, current_locus$lower_pos_lod_drop) * 1e6
    locus_end_bp <- max(current_locus$upper_pos_lod_drop, current_locus$lower_pos_lod_drop) * 1e6
    
    # Check which genes in our main table are within this specific locus
    main_genes_table[, (locus_col_name) := "No"]  # Default all to FALSE
    genes_in_this_locus <- main_genes_table[chr == current_locus$chr & 
                                              start_bp < locus_end_bp & 
                                              end_bp > locus_start_bp, 
                                            `Mouse Ensembl ID`]
    
    # Set to TRUE for genes in this locus
    main_genes_table[`Mouse Ensembl ID` %in% genes_in_this_locus, (locus_col_name) := "Yes"]
  }
  
  # Add the locus columns to our column ordering
  locus_cols <- grep("^In_", names(main_genes_table), value = TRUE)
  remaining_cols <- setdiff(names(main_genes_table), c(initial_cols, locus_cols))
  
  # Create final column order: initial columns, locus presence columns, then remaining columns
  col_order <- c(initial_cols, locus_cols, remaining_cols)
  main_genes_table <- main_genes_table[, ..col_order]
  
  # --- Create Excel Workbook ---
  wb <- createWorkbook()
  
  # Add the main summary sheet with all genes and locus presence
  addWorksheet(wb, "AllGenesInCluster")
  writeData(wb, "AllGenesInCluster", main_genes_table)
  
  # Style the locus presence columns as TRUE/FALSE
  locus_col_indices <- which(names(main_genes_table) %in% locus_cols)
  for (col_idx in locus_col_indices) {
    conditionalFormatting(wb, "AllGenesInCluster", 
                          cols = col_idx, 
                          rows = 1:(nrow(main_genes_table) + 1), # +1 for header
                          rule = "==Yes", 
                          style = createStyle(bgFill = "#90EE90")) # Light green for TRUE
  }
  
  # Add sheets for individual genes with phenotype/disease data
  for (idx in 1:nrow(main_genes_table)) {
    current_gene_symbol <- main_genes_table$`Mouse Gene Symbol`[idx]
    current_mouse_ensembl_id <- main_genes_table$`Mouse Ensembl ID`[idx]
    
    # Check for and add Mouse Phenotypes sheet
    gene_mouse_pheno_data <- full_mouse_pheno_data[mouse_gene_symbol == current_gene_symbol, 
                                                   .(mouse_gene_symbol, 
                                                     OntologyTermID = `OntologyAnnotation.ontologyTerm.identifier`, 
                                                     OntologyTermName = `OntologyAnnotation.ontologyTerm.name`,
                                                     PubMedID = `OntologyAnnotation.evidence.publications.pubMedId`,
                                                     Description = `OntologyAnnotation.evidence.comments.description`)]
    
    if (nrow(gene_mouse_pheno_data) > 0) {
      sheet_name_mouse <- substr(paste0(current_gene_symbol, "_MousePheno"), 1, 31) # Excel sheet name limit
      addWorksheet(wb, sheet_name_mouse)
      writeData(wb, sheet_name_mouse, gene_mouse_pheno_data)
    }
    
    # Check for and add Human Disease Associations sheet
    # First, get human ortholog info
    human_ortho_info <- orthology_data[mouse_ensembl_id == current_mouse_ensembl_id, ]
    
    if (nrow(human_ortho_info) > 0 && human_ortho_info$human_ensembl_id[1] != "") {
      current_human_ensembl_id <- human_ortho_info$human_ensembl_id[1]
      current_human_symbol <- human_ortho_info$human_gene_symbol[1]
      
      gene_human_disease_data <- full_human_assoc_data[human_ensembl_id == current_human_ensembl_id, 
                                                       .(human_ensembl_id, 
                                                         human_gene_symbol = symbol, 
                                                         disease_id, 
                                                         disease_name, 
                                                         association_score)]
      
      if (nrow(gene_human_disease_data) > 0) {
        sheet_name_human <- substr(paste0(current_gene_symbol, "_HumanDisease"), 1, 31)
        addWorksheet(wb, sheet_name_human)
        writeData(wb, sheet_name_human, gene_human_disease_data)
      }
    }
  }
  
  # --- Save Workbook ---
  # Use the first locus info for naming
  first_locus <- loci_info_df[1, ]
  excel_file_name <- file.path(output_dir, paste0("gene_info_cluster_chr", first_locus$chr, "_", 
                                                  round(min(loci_info_df$upper_pos_lod_drop), 0), "-", 
                                                  round(max(loci_info_df$lower_pos_lod_drop), 0), "Mb.xlsx"))
  
  message("Generating gene information Excel file: ", excel_file_name)
  saveWorkbook(wb, excel_file_name, overwrite = TRUE)
  
  return(main_genes_table)
} 


# Function to generate Founder Mutation Table (for SNPs in locus)
generate_founder_mutation_table <- function(founder_snps_in_locus_dt, # data.table from get_founder_snps_in_region
                                            locus_info,
                                            output_dir) {
  if (nrow(founder_snps_in_locus_dt) == 0) {
    message("No founder SNPs found in the specified locus for table.")
    return(NULL)
  }
  
  founder_snp_table <- founder_snps_in_locus_dt[, .(
    SNP_ID = SNP,
    rsID,
    CHR,
    Position_mm39 = mm39,
    MajorAllele = Major,
    MinorAllele = Minor,
    MAF,
    Strains_with_Minor_Allele = SNPStrains
    # Add columns for consequence, overlapping gene if you implement that
  )]
  
  table_file_name <- file.path(output_dir, paste0("founder_snp_table_chr", locus_info$chr, "_", round(locus_info$peak_pos,2), "Mb.csv"))
  message("Generating founder SNP table: ", table_file_name)
  fwrite(founder_snp_table, table_file_name)
  return(table_file_name)
}



# Create the summary text with more detailed contextual information
summary_text <- paste( # this prints s
  "# CARDIAC GENE SUMMARY\n",
  "Summary for Locus:", locus_packet_name,
  "\nAssociated Trait(s):", 
  "\n    Ctrl: ", paste(loci |> filter(drug == "Ctrl") |>  pull(trait), collapse= ", "),
  "\n    Iso:", paste(loci |> filter(drug == "Iso") |>  pull(trait), collapse= ", "),    
  "\nGenomic Region (chr", current_locus_info$chr, "): ", min(loci$upper_pos_lod_drop), "-", max(loci$lower_pos_lod_drop), " (mm39)",
  "\n\n# Summary of known assocations of genes with this loci",
  "\nTotal cardiac-related genes in region: ", total_count,
  "\n  - Human cardiac traits only: ", human_only_count,
  "\n  - Mouse cardiac phenotypes only: ", mouse_only_count,
  "\n  - Both human and mouse: ", both_count,
  "\n\n# DETAILED GENE ANNOTATIONS",
  "\n\n## 1. Genes known to be associated with human cardiac traits only (", human_only_count, "):",
  "\n        ", format_human_only(human_only_genes),
  "\n\n## 2. Genes known to be associated with mouse cardiac phenotypes only (", mouse_only_count, "):",
  "\n        ", format_mouse_only(mouse_only_genes),
  "\n\n## 3. Genes associated with both human and mouse cardiac traits (", both_count, "):",
  "\n        ", format_both(human_mouse_genes),
  "\n\n# NOTES",
  "\n- All 'Associated Traits (Drug)' in the original data are our own miQTL mappings",
  "\n- See cardiac_genes_summary.csv for a structured version of this data"
)

# Write to README_summary.txt
writeLines(summary_text, file.path(current_output_dir, "README_summary.txt"))
# Create a more structured CSV output with additional useful metadata
cardiac_genes_csv <- rbind(
  if(nrow(human_only_genes) > 0) {
    data.frame(
      Gene_Symbol = human_only_genes$`Mouse Gene Symbol`,
      Ensembl_ID = human_only_genes$`Mouse Ensembl ID`,
      Category = "Human Only",
      Human_Cardiac_Traits = sapply(human_only_genes$`Human Disease Summary`, extract_human_cardiac_traits),
      Mouse_Cardiac_Phenotypes = NA,
      Chromosome = human_only_genes$chr,
      Start_Position = human_only_genes$start_bp,
      End_Position = human_only_genes$end_bp,
      NRVM_CPM= human_only_genes$`Avg NRVM CPM (Ctrl)`,
      stringsAsFactors = FALSE
    )
  },
  if(nrow(mouse_only_genes) > 0) {
    data.frame(
      Gene_Symbol = mouse_only_genes$`Mouse Gene Symbol`,
      Ensembl_ID = mouse_only_genes$`Mouse Ensembl ID`,
      Category = "Mouse Only",
      Human_Cardiac_Traits = NA,
      Mouse_Cardiac_Phenotypes = sapply(mouse_only_genes$`Mouse Phenotype Summary (MGI)`, extract_mouse_cardiac_traits),
      Chromosome = mouse_only_genes$chr,
      Start_Position = mouse_only_genes$start_bp,
      End_Position = mouse_only_genes$end_bp,
      NRVM_CPM = mouse_only_genes$`Avg NRVM CPM (Ctrl)`,
      stringsAsFactors = FALSE
    )
  },
  if(nrow(human_mouse_genes) > 0) {
    data.frame(
      Gene_Symbol = human_mouse_genes$`Mouse Gene Symbol`,
      Ensembl_ID = human_mouse_genes$`Mouse Ensembl ID`,
      Category = "Both",
      Human_Cardiac_Traits = sapply(human_mouse_genes$`Human Disease Summary`, extract_human_cardiac_traits),
      Mouse_Cardiac_Phenotypes = sapply(human_mouse_genes$`Mouse Phenotype Summary (MGI)`, extract_mouse_cardiac_traits),
      Chromosome = human_mouse_genes$chr,
      Start_Position = human_mouse_genes$start_bp,
      End_Position = human_mouse_genes$end_bp,
      NRVM_CPM = human_mouse_genes$`Avg NRVM CPM (Ctrl)`,
      stringsAsFactors = FALSE
    )
  }
)

# Also create an Excel file for better visualization and filtering capabilities
if(require(openxlsx)) {
  wb <- createWorkbook()
  
  # Add a summary sheet
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", data.frame(
    Metric = c(
      "Locus", 
      "Associated Trait", 
      "Genomic Region", 
      "Peak Position (miQTL)", 
      "Max LOD (miQTL)",
      "Total Cardiac Genes",
      "Human Only",
      "Mouse Only",
      "Both Human and Mouse"
    ),
    Value = c(
      locus_packet_name,
      current_locus_info$trait,
      paste0("chr", current_locus_info$chr, ": ", locus_start_bp, "-", locus_end_bp, " (mm39)"),
      paste0(round(current_locus_info$peak_pos, digits = 2), " Mb"),
      round(current_locus_info$max_lod, digits = 1),
      total_count,
      human_only_count,
      mouse_only_count,
      both_count
    )
  ))
  
  # Add a gene details sheet
  addWorksheet(wb, "Cardiac Genes")
  writeData(wb, "Cardiac Genes", cardiac_genes_csv)
  
  # Style the workbook
  setColWidths(wb, "Summary", cols = 1:2, widths = c(30, 50))
  setColWidths(wb, "Cardiac Genes", cols = 1:9, widths = c(15, 25, 15, 50, 50, 10, 15, 15, 15))
  
  # Save the workbook
  saveWorkbook(wb, file.path(current_output_dir, "cardiac_genes_summary.xlsx"), overwrite = TRUE)
}