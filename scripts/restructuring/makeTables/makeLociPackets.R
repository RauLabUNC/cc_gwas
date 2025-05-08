# --- 0. Load Libraries ---
# Core libraries
library(tidyverse)
library(data.table)
library(plotgardener) # For locus plots
library(TxDb.Mmusculus.UCSC.mm39.knownGene) # For gene models if needed by plotgardener
library(org.Mm.eg.db) # For gene ID mapping if needed
library(openxlsx) # For writing Excel files

# --- 1. Load All Data --- ####
# (Using the paths and loading logic you provided)

message("Loading core data tables...")

# === From multTrait...Exp.R ===
base_relational <- "data/processed/joinLoci/relational_tables"
paths_multitrait <- list(
  genes_mouse   = fs::path(base_relational, "genes_mouse.csv"),
  orthology     = fs::path(base_relational, "orthology.csv"),
  associations  = fs::path(base_relational, "associations.csv"),
  trait_loci    = fs::path(base_relational, "traitLoci.csv"),
  exp_loci_cis  = fs::path(base_relational, "expLoci_cisflag.csv"),
  mouse_pheno   = fs::path(base_relational, "mouseGenePhenotypes.csv")
  # nrvm_counts and nrvm_meta are also in your original script, load if direct NRVM expression plots are needed
)

genes_mouse   <- fread(paths_multitrait$genes_mouse)
ortho_mouse2h <- fread(paths_multitrait$orthology)
trait_loci    <- fread(paths_multitrait$trait_loci) # Contains miQTL trait loci
exp_cis       <- fread(paths_multitrait$exp_loci_cis)
associations  <- fread(paths_multitrait$associations) # Human disease associations
mouse_pheno   <- fread(paths_multitrait$mouse_pheno)

# === Raw and normalized phenotypes for your study ===
# Correcting file names based on your loading script (Ctrl data from Iso file and vice-versa)
# Assuming the file names in your script were swapped:
pheno_box_ctrl <- fread("data/processed/phenotypes/boxCoxTest/boxcox_individual_Ctrl.csv") # Should contain Drug == "Ctrl"
pheno_box_iso  <- fread("data/processed/phenotypes/boxCoxTest/boxcox_individual_Iso.csv")  # Should contain Drug == "Iso"
# It's good practice to combine these and filter by Drug column later
all_pheno_data <- rbind(pheno_box_ctrl, pheno_box_iso, fill = TRUE) # fill=TRUE if columns might differ slightly

# === Data used in plotLoci ===
sig_regions    <- read_csv("data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv") %>%
  mutate(
    across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
    chr = as.character(chr)
  ) %>%
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug)

scan_data      <- readRDS("results/sig_regions/scan_data.rds") # miQTL scan results
threshold_data <- readRDS("results/sig_regions/threshold_data.rds")

pyResults <- list() # Changed from c() to list() for proper assignment
pyResults[["Ctrl"]] <- fread("data/processed/joinLoci/trait_qtl/PyLMM/Ctrl_pvals.csv")
pyResults[["Iso"]]  <- fread("data/processed/joinLoci/trait_qtl/PyLMM/Iso_pvals.csv")

# === Merged table output from multTrait_cis-eQTL_nrvmExp.R ===
# This table is very useful as it pre-joins a lot of relevant gene info
merged_gene_info <- fread("results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv")

# === CC Founder Variant Information ===
muNoSplice <- fread("data/processed/joinLoci/relational_tables/ccVariants/Gene_Mutations_CC_WT_NoSplice_OnlyDeleterious_3_31_25.csv")
muDriving  <- fread("data/processed/joinLoci/relational_tables/ccVariants/Driving_Mutations_by_SNP_Filtered.csv")


# === get the genome assembly

# Genome assembly
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

message("All data loaded.")

# --- 2. Define Helper Functions --- ####

# Helper to get genes within a given genomic region
get_genes_in_region <- function(target_chr, target_start, target_end, gene_data = genes_mouse) {
  gene_data[chr == target_chr & start_bp < target_end & end_bp > target_start, ]
}

get_high_evidence_genes_in_region <-function(target_chr, target_start, target_end, gene_data = merged_gene_info){
  merged_gene_info |> filter(mouse_gene_symbol %in% genes_in_current_locus$mouse_gene_symbol) |> pull(mouse_gene_symbol)
}
# Helper to get founder SNPs in a given genomic region
get_founder_snps_in_region <- function(target_chr, target_start, target_end, founder_snp_data = muDriving) {
  founder_snp_data[CHR == target_chr & mm39 < target_end & mm39 > target_start, ]
}

# Helper to get gene-specific founder mutations (from muNoSplice)
get_gene_founder_mutations <- function(gene_symbol, mutation_data = muNoSplice) {
  mutation_data[Gene == gene_symbol, ]
}


# --- 3. Define Plotting Functions --- ####

# Function to generate Locus Zoom Plot (simplified from your plotLoci_pyLMMvsMiQTL.R)
# This would be a more complex function, reusing logic from your existing script
generate_locus_zoom_plot <- function(locus_info, # A row from sig_regions or a custom list
                                     output_dir,
                                     genes_in_locus,
                                     top_genes_in_locus,
                                     founder_snps_in_locus) {
  # locus_info should contain: chr, peak_pos, lower_pos_lod_drop, upper_pos_lod_drop, trait, drug
  #locus_info <- current_locus_info # A row from sig_regions or a custom list
  #output_dir <- current_output_dir
  #genes_in_locus <- genes_in_current_locus
  #top_genes_in_locus <- top_genes_in_locus
  #founder_snps_in_locus <- founder_snps_current_locus
  
  
  # Define plot parameters (chromosome, start, end)
  current_chr <- locus_info$chr
  # Ensure positions are in base pairs for plotgardener
  # Your sig_regions peak_pos etc. are in Mb, so multiply by 1e6
  # Padded region for plotting:
  plot_start_bp <- floor(locus_info$upper_pos_lod_drop * 1e6) - 5e5 # Add padding
  plot_end_bp   <- ceiling(locus_info$lower_pos_lod_drop * 1e6) + 5e5 # Add padding
  
  bounds_bp <- c(locus_info$upper_pos_lod_drop  * 1e6, locus_info$lower_pos_lod_drop * 1e6)
  # Ensure start is not negative
  plot_start_bp <- max(0, plot_start_bp)
  
  
  plot_file_name <- file.path(output_dir,
                              paste0("locus_zoom_", locus_info$trait, "_", locus_info$drug,
                                     "_chr", current_chr, "_", round(locus_info$peak_pos, 2), "Mb.pdf"))
  
  message("Generating locus zoom plot: ", plot_file_name)
  
  # --- plotgardener setup ---
  pdf(plot_file_name, width = 9, height = 5.5)
  pageCreate(width = 9, height = 5.5, default.units = "inches", showGuides = FALSE)
  
  params_genome <- pgParams(
    assembly   = mm39, 
    chrom      = paste0("chr", current_chr),
    chromstart = plot_start_bp,
    chromend   = plot_end_bp
  )
  
  # --- Plot miQTL LOD scores ---
  # (Adapted from your plotLoci_pyLMMvsMiQTL.R)
  # Requires scan_data for the current trait and drug
  current_scan_data <- scan_data[[locus_info$trait]][[locus_info$drug]]
  miqtl_df_for_plot <- tibble(
    marker = names(current_scan_data$LOD),
    chr    = as.character(current_scan_data$chr),
    pos    = current_scan_data$pos$Mb * 1e6, # Ensure pos is in bp
    lod    = current_scan_data$LOD
  ) %>%
    filter(chr == current_chr, !is.na(pos), !is.na(lod)) %>% # Filter for current chromosome
    transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod))
  
  # Determine y-axis limits for miQTL plot
  # Placeholder for threshold - you'd get this from your threshold_data
  miqtl_threshold_val <- threshold_data[[locus_info$trait]][[locus_info$drug]] # Example threshold
  miqtl_ylim <- c(0, max(c(-log10(miqtl_df_for_plot$p), miqtl_threshold_val, 5), na.rm = TRUE) + 1)
  
  miqtl_plot <- plotManhattan(
    data = miqtl_df_for_plot, params = params_genome,
    range = miqtl_ylim, trans = "-log10", sigVal = 10^(-miqtl_threshold_val),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y ,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), xfield = "pos", yfield = "p",
    fill = "#a6cee3", sigCol = "#1f78b4", sigLine = TRUE, baseline = TRUE,
    default.units = "inches"
  )
  
  annoYaxis(plot = miqtl_plot, at = pretty(miqtl_ylim),
            axisLine     = TRUE,
            fontsize     = 8,
            main         = FALSE)
  plotText(label = "LOD (miQTL)",  
           x       = 2 * PLOT_PARAMS$x + 0.1,
           y       = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height / 2,
           rot     = 270,
           fontsize= 8,
           just    = "center",
           default.units = "in")
  # --- Plot pyLMM p-values ---
  # (Adapted from your plotLoci_pyLMMvsMiQTL.R)
  current_pylmm_data <- pyResults[[locus_info$drug]]
  pylmm_df_for_plot <- current_pylmm_data %>%
    filter(Chr == current_chr) %>%
    transmute(
      chrom = paste0("chr", Chr),
      pos   = Pos_mm39, # Already in bp
      p = .data[[locus_info$trait]] # p-values are in trait columns
    ) %>% filter(!is.na(p) & p > 0) # Ensure p is valid for -log10
  
  pylmm_ylim <- c(0, ceiling(max(-log10(pylmm_df_for_plot$p), 5, na.rm = TRUE)) +1)
  
  pylmm_plot <- plotManhattan(
    data = pylmm_df_for_plot, params = params_genome,
    range = pylmm_ylim, trans = "-log10",
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height + 0.2,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), xfield = "pos", yfield = "p", # yfield is p_val
    fill = "black", baseline = TRUE, sigLine = FALSE, # No threshold line for pyLMM by default in your script
    default.units = "inches"
  )
  annoYaxis(plot = pylmm_plot, at = pretty(pylmm_ylim),
            axisLine     = TRUE,
            fontsize     = 8,
            main         = FALSE)
  plotText(label = "-log10(p) (pyLMM)",
           x       = 2 * PLOT_PARAMS$x + 0.1,
           y       = PLOT_PARAMS$plot_y + 1.5*PLOT_PARAMS$plot_height + 0.2,
           rot     = 270,
           fontsize= 8,
           just    = "center",
           default.units = "in")
  # Highlight significant region on main_plot
  annoHighlight(
    plot         = miqtl_plot,
    chrom        = paste0("chr", current_chr),
    chromstart   = floor(min(bounds_bp)),
    chromend     = ceiling(max(bounds_bp)),
    fill         = "#fb9a99",
    y            = PLOT_PARAMS$plot_y,
    height       = PLOT_PARAMS$plot_height,
    just         = c("left", "top"),
    default.units= "inches",
    alpha        = 0.2,
    params       = params_genome
  )
  
  # --- Plot Genes ---
  # This requires `genes_in_locus` to be prepared
  # For plotGenes, data needs to be in a specific format or use TxDb
  # Simplification: using TxDb directly if genes_in_locus is not pre-formatted for plotGenes
  genes_y_pos <- 0.5 + 1.5 + 0.2 + 1.5 + 0.2 # Position below pyLMM
  
  gene_plot <- plotGenes(
    params = params_genome,
    x = 4.25, y = genes_y_pos, width = 8, height = 1, # Adjust height as needed
    just = c("center", "top"), default.units = "inches",
    geneOrder = "score", # Example, or remove for default ordering
    fontsize = 7,
    geneHighlights = top_genes_in_locus,
    geneBackground = "#fdbb84"
  )
  
  # --- Plot Founder Variants (if data available) ---
  
  
  # --- Add Genome Label ---
  annoGenomeLabel(
    plot = gene_plot, # Attach to one of the plots, e.g., gene_plot
    params = params_genome,
    x = 4.25, y = genes_y_pos + 1, # Below last track
    scale = "Mb", fontsize = 10,
    just = c("center", "top"), default.units = "inches"
  )
  
  # --- Add Title ---
  plotText(
    label = paste("QTL:", locus_info$trait, "(",locus_info$drug, ") - Chr",
                  current_chr, "Peak:", round(locus_info$peak_pos, 2), "Mb"),
    x = 4.5, y = 0.1, just = c("center", "top"),
    fontface = "bold", fontsize = 12, default.units = "inches"
  )
  
  dev.off()
  return(plot_file_name)
}

# Function to generate Phenotype Distribution Plot
generate_phenotype_dist_plot <- function(trait_name, drug_condition,
                                         pheno_data = all_pheno_data,
                                         output_dir) {
  
  # Filter for the specific drug condition
  current_pheno_data <- pheno_data[Drug == drug_condition, ]
  
  # Check if trait_name exists in the data
  if (!trait_name %in% names(current_pheno_data)) {
    warning(paste("Trait", trait_name, "not found in phenotype data for", drug_condition))
    return(NULL)
  }
  
  plot_file_name <- file.path(output_dir, paste0("pheno_dist_", trait_name, "_", drug_condition, ".png"))
  message("Generating phenotype distribution plot: ", plot_file_name)
  
  # Create the plot (example: density plot)
  # This assumes your pheno_box_... files contain the raw (or Box-Cox transformed) values
  p <- ggplot(current_pheno_data, aes_string(x = trait_name)) +
    geom_density(fill = "skyblue", alpha = 0.7) +
    geom_histogram(aes(y = ..density..), binwidth = NULL, fill="grey", alpha=0.5, color="black") + # Added histogram
    labs(title = paste("Distribution of", trait_name, "for", drug_condition),
         x = trait_name, y = "Density") +
    theme_minimal()
  
  ggsave(plot_file_name, plot = p, width = 6, height = 4)
  return(plot_file_name)
}

# Function to generate Gene Expression Plot (e.g., for a specific gene)
# This would require cardiac tissue expression data (not just NRVM averages if plotting by group)
# For now, let's assume you might have a table `cardiac_expression_data` with columns:
# mouse_gene_symbol, SampleID, expression_value, Drug (Ctrl/Iso)
generate_gene_expression_plot <- function(gene_symbol,
                                          # cardiac_exp_data, # This data needs to be loaded/defined
                                          output_dir) {
  # Placeholder: This function needs actual expression data to work.
  # if (missing(cardiac_exp_data) || !gene_symbol %in% cardiac_exp_data$mouse_gene_symbol) {
  #   message("Cardiac expression data not available or gene not found for: ", gene_symbol)
  #   return(NULL)
  # }
  #
  # gene_specific_exp <- cardiac_exp_data[mouse_gene_symbol == gene_symbol, ]
  #
  # plot_file_name <- file.path(output_dir, paste0("gene_exp_", gene_symbol, ".png"))
  # message("Generating gene expression plot: ", plot_file_name)
  #
  # p <- ggplot(gene_specific_exp, aes(x = Drug, y = expression_value, fill = Drug)) +
  #   geom_boxplot(outlier.shape = NA) +
  #   geom_jitter(width = 0.1, alpha = 0.5) +
  #   labs(title = paste("Expression of", gene_symbol, "in Cardiac Tissue"),
  #        x = "Condition", y = "Normalized Expression") +
  #   theme_minimal()
  #
  # ggsave(plot_file_name, plot = p, width = 5, height = 5)
  # return(plot_file_name)
  message("Placeholder for gene expression plot for: ", gene_symbol)
  return(NULL)
}


# --- 4. Define Table Generation Functions --- ####

# Function to generate Gene Information Excel File with multiple sheets
generate_gene_info_excel <- function(genes_in_locus_dt, 
                                     locus_info, 
                                     output_dir,
                                     full_mouse_pheno_data = mouse_pheno,
                                     full_human_assoc_data = associations,
                                     orthology_data = ortho_mouse2h) {
  
  if (nrow(genes_in_locus_dt) == 0) {
    message("No genes found in the specified locus for Excel generation.")
    return(NULL)
  }
  
  # --- Create the main summary table (similar to before) ---
  ensembl_ids_in_locus <- genes_in_locus_dt$mouse_ensembl_id
  locus_gene_summary <- merged_gene_info[mouse_ensembl_id %in% ensembl_ids_in_locus, ]
  
  founder_muts_gene_level <- muNoSplice[Gene %in% locus_gene_summary$mouse_gene_symbol, 
                                        .(mouse_gene_symbol = Gene, founder_gene_mutations = Mutations)]
  
  if (nrow(founder_muts_gene_level) > 0) {
    locus_gene_summary <- merge(locus_gene_summary, founder_muts_gene_level, by = "mouse_gene_symbol", all.x = TRUE)
  } else {
    locus_gene_summary[, founder_gene_mutations := NA_character_]
  }
  
  main_summary_table <- locus_gene_summary[, .(
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
  
  gene_coords <- genes_mouse[mouse_ensembl_id %in% main_summary_table$`Mouse Ensembl ID`, 
                             .(mouse_ensembl_id, chr, start_bp, end_bp)]
  setnames(gene_coords, "mouse_ensembl_id", "Mouse Ensembl ID")
  main_summary_table <- merge(main_summary_table, gene_coords, by = "Mouse Ensembl ID", all.x = TRUE)
  
  col_order <- c("Mouse Gene Symbol", "Mouse Ensembl ID", "chr", "start_bp", "end_bp", 
                 setdiff(names(main_summary_table), c("Mouse Gene Symbol", "Mouse Ensembl ID", "chr", "start_bp", "end_bp")))
  main_summary_table <- main_summary_table[, ..col_order]
  
  # --- Create Excel Workbook ---
  wb <- createWorkbook()
  addWorksheet(wb, "LocusGeneSummary")
  writeData(wb, "LocusGeneSummary", main_summary_table)
  
  # --- Add sheets for individual genes with phenotype/disease data ---
  for (idx in 1:nrow(main_summary_table)) {
    current_gene_symbol <- main_summary_table$`Mouse Gene Symbol`[idx]
    current_mouse_ensembl_id <- main_summary_table$`Mouse Ensembl ID`[idx]
    
    # Check for and add Mouse Phenotypes sheet
    gene_mouse_pheno_data <- full_mouse_pheno_data[mouse_gene_symbol == current_gene_symbol, 
                                                   .(mouse_gene_symbol, 
                                                     OntologyTermID = `OntologyAnnotation.ontologyTerm.identifier`, # Renaming for clarity
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
                                                         association_score)] # Add all relevant columns from 'associations'
      
      if (nrow(gene_human_disease_data) > 0) {
        sheet_name_human <- substr(paste0(current_gene_symbol, "_HumanDisease"), 1, 31)
        addWorksheet(wb, sheet_name_human)
        writeData(wb, sheet_name_human, gene_human_disease_data)
      }
    }
  }
  
  # --- Save Workbook ---
  excel_file_name <- file.path(output_dir, paste0("gene_info_details_chr", locus_info$chr, "_", round(locus_info$peak_pos,2), "Mb.xlsx"))
  message("Generating gene information Excel file: ", excel_file_name)
  saveWorkbook(wb, excel_file_name, overwrite = TRUE)
  
  return(excel_file_name)
}

# Function to generate Founder Mutation Table (for SNPs in locus)
generate_founder_mutation_table <- function(founder_snps_in_locus_dt, # data.table from get_founder_snps_in_region
                                            locus_info,
                                            output_dir) {
  if (nrow(founder_snps_in_locus_dt) == 0) {
    message("No founder SNPs found in the specified locus for table.")
    return(NULL)
  }
  
  # Select relevant columns from muDriving
  # You might want to link these SNPs to genes if they fall within gene boundaries
  # For now, just a table of SNPs in the locus
  
  # Identify genes overlapping these SNPs
  # This is a bit more involved: iterate through SNPs and find overlapping genes
  # For simplicity, we'll just output the SNP table.
  
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


# --- 5. Main Loop to Generate Packets --- ####

# Define which loci to process. For example, iterate through `sig_regions`
# Or define a custom list of loci or genes.
# For testing, let's pick the first few significant regions
loci_to_process <- sig_regions[11, ] # Process first 2 significant regions for example

# Base output directory for all packets
base_output_dir <- "results/qtl_packets" 
if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

for (i in 1:nrow(loci_to_process)) {
  current_locus_info <- loci_to_process[i, ]
  
  # Create a specific directory for this locus packet
  locus_packet_name <- paste0("locus_chr", current_locus_info$chr,
                              "_peak", round(current_locus_info$peak_pos, 2), "Mb",
                              "_", current_locus_info$trait, "_", current_locus_info$drug)
  current_output_dir <- file.path(base_output_dir, locus_packet_name)
  if (!dir.exists(current_output_dir)) dir.create(current_output_dir, recursive = TRUE)
  
  message(paste0("\n--- Processing Locus: ", locus_packet_name, " ---"))
  
  # --- A. Get genes and founder SNPs in the current locus region ---
  # Define locus boundaries (e.g., from LOD drop or fixed window around peak)
  # Using LOD drop from sig_regions, converting Mb to bp
  locus_start_bp <- ceiling(current_locus_info$upper_pos_lod_drop * 1e6)
  locus_end_bp   <- floor(current_locus_info$lower_pos_lod_drop * 1e6)
  
  genes_in_current_locus <- get_genes_in_region(current_locus_info$chr,
                                                locus_start_bp,
                                                locus_end_bp)
  
  founder_snps_current_locus <- get_founder_snps_in_region(current_locus_info$chr,
                                                           locus_start_bp,
                                                           locus_end_bp)
  
  top_genes_in_locus <- data.frame(gene = get_high_evidence_genes_in_region(),
                                   color = "#e34a33")
  
  
  # --- B. Generate Locus Zoom Plot ---
  # Constants
  COLORS <- c(Sig = "#1f78b4", nonSig = "#a6cee3")
  PLOT_DIMS <- list(page_width = 9, page_height = 6.5, res = 300)
  PLOT_PARAMS <- list(x = 4.25, plot_width = 8, plot_height = 1.5, plot_y = 0.5)
  GENE_DIMS <- list(height = 2, y_offset = 0.5, label_offset = 0.1)
  MIN_YLIM <- 5
  
  # --- B. Generate Locus Zoom Plot ---
  locus_plot_file <- generate_locus_zoom_plot(current_locus_info,
                                              current_output_dir,
                                              genes_in_current_locus, 
                                              top_genes_in_locus, # Pass this for potential use in plotGenes
                                              founder_snps_current_locus) # Pass for founder variant track
  
  # --- C. Generate Phenotype Distribution Plot ---
  pheno_dist_plot_file <- generate_phenotype_dist_plot(current_locus_info$trait,
                                                       current_locus_info$drug,
                                                       pheno_data = all_pheno_data,
                                                       output_dir = current_output_dir)
  
  # --- D. Generate Gene Information Table ---
  gene_info_table_file <- generate_gene_info_table(genes_in_current_locus,
                                                   current_locus_info,
                                                   current_output_dir)
  
  # UPDATED CALL to generate Excel file
  gene_info_excel_file <- generate_gene_info_excel(genes_in_current_locus,
                                                   current_locus_info,
                                                   current_output_dir,
                                                   full_mouse_pheno_data = mouse_pheno,
                                                   full_human_assoc_data = associations,
                                                   orthology_data = ortho_mouse2h)
  
  # --- E. Generate Founder Mutation Table for SNPs in Locus ---
  founder_snp_table_file <- generate_founder_mutation_table(founder_snps_current_locus,
                                                            current_locus_info,
                                                            current_output_dir)
  
  # --- F. Generate Gene Expression Plots (for top N genes or specific candidates) ---
  # Example: Plot expression for the first gene in the locus, if any
  if (nrow(genes_in_current_locus) > 0) {
    candidate_gene_symbol <- genes_in_current_locus$mouse_gene_symbol[1]
    # gene_exp_plot_file <- generate_gene_expression_plot(candidate_gene_symbol,
    #                                                     # cardiac_expression_data = your_cardiac_exp_data, # You need to load this
    #                                                     output_dir = current_output_dir)
  }
  
  # --- G. (Manual Step Idea) Create a Summary README.txt for the packet ---
  summary_text <- paste(
    "Summary for Locus:", locus_packet_name,
    "\nAssociated Trait:", current_locus_info$trait,
    "\nCondition:", current_locus_info$drug,
    "\nGenomic Region (chr", current_locus_info$chr, "): ", locus_start_bp, "-", locus_end_bp, " (mm39)",
    "\nPeak Position (miQTL):", current_locus_info$peak_pos, "Mb",
    "\nMax LOD (miQTL):", current_locus_info$max_lod,
    "\n\nFiles in this packet:",
    if (!is.null(locus_plot_file)) paste0("\n- Locus Zoom Plot: ", basename(locus_plot_file)),
    if (!is.null(pheno_dist_plot_file)) paste0("\n- Phenotype Distribution Plot: ", basename(pheno_dist_plot_file)),
    if (!is.null(gene_info_table_file)) paste0("\n- Gene Information Table: ", basename(gene_info_table_file)),
    if (!is.null(founder_snp_table_file)) paste0("\n- Founder SNP Table: ", basename(founder_snp_table_file))
    # if (!is.null(gene_exp_plot_file)) paste0("\n- Example Gene Expression Plot: ", basename(gene_exp_plot_file))
  )
  writeLines(summary_text, file.path(current_output_dir, "README_summary.txt"))
  
  message(paste0("--- Finished processing packet for: ", locus_packet_name, " ---"))
}

message("\nAll selected loci processed. Packets generated in: ", base_output_dir)

