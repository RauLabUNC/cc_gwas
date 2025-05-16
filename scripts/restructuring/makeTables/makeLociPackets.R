# --- 0. Load Libraries ---
# Core libraries
library(tidyverse)
library(data.table)
library(plotgardener) # For locus plots
library(TxDb.Mmusculus.UCSC.mm39.knownGene) 
library(org.Mm.eg.db) # For gene ID mapping
library(openxlsx) # For writing Excel files
library(miqtl)
library(igraph)          # install.packages("igraph") if you don’t have it
library(GenomicRanges)
library(RColorBrewer)

# --- 1. Load All Data --- ####
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
generate_locus_zoom_plot <- function(locus_info, # A row from sig_regions or a custom list
                                     output_dir,
                                     genes_in_locus,
                                     top_genes_in_locus,
                                     founder_snps_in_locus) {
  # locus_info should contain: chr, peak_pos, lower_pos_lod_drop, upper_pos_lod_drop, trait, drug
  locus_info <- current_locus_info # A row from sig_regions or a custom list
  output_dir <- current_output_dir
  genes_in_locus <- genes_in_current_locus
  top_genes_in_locus <- top_genes_in_locus
  founder_snps_in_locus <- founder_snps_current_locus
  
  # Define plot parameters (chromosome, start, end)
  current_chr <- locus_info$chr
  
  # Ensure positions are in base pairs for plotgardener
  plot_start_bp <- floor(locus_info$upper_pos_lod_drop * 1e6) - 1e6 # Add padding
  plot_end_bp   <- ceiling(locus_info$lower_pos_lod_drop * 1e6) + 1e6 # Add padding
  bounds_bp <- c(locus_info$upper_pos_lod_drop  * 1e6, locus_info$lower_pos_lod_drop * 1e6)
  
  # Ensure start is not negative
  plot_start_bp <- max(0, plot_start_bp)
  
  dir_path <- file.path(output_dir, "zoomPlots")
  plot_file_name <- file.path(dir_path,
                              paste0("locus_zoom_", locus_info$trait, "_", locus_info$drug,
                                     "_chr", current_chr, "_", round(locus_info$peak_pos, 2), "Mb.pdf"))
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  message("Generating locus zoom plot: ", plot_file_name)
  
  # --- plotgardener setup ---
  pdf(plot_file_name, width = 10.5, height = 6.5)
  pageCreate(width = 10.5, height = 6.5, default.units = "inches", showGuides = FALSE)

  params_genome <- pgParams(
    assembly   = mm39, 
    chrom      = paste0("chr", current_chr),
    chromstart = plot_start_bp,
    chromend   = plot_end_bp
  )
  
  # --- Plot miQTL LOD scores ---
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
  
  # Plot allele effects
  
  pos_in_range_logical <- current_scan_data$pos$Mb * 10^6 >= plot_start_bp & current_scan_data$pos$Mb * 10^6 <= plot_end_bp
  chr_logical <- current_scan_data$chr == current_locus_info$chr
  marker_positions_bp <- current_scan_data$pos$Mb[pos_in_range_logical & chr_logical] * 10^6
  
  markers <- data.frame(marker = current_scan_data$loci[pos_in_range_logical & chr_logical],
                        start  = marker_positions_bp) |> 
    arrange(start)
  
  allele_effects_matrix <- current_scan_data$allele.effects[,markers$marker, drop = FALSE] # drop=FALSE to handle single marker case
  allele_effects_transposed <- t(current_scan_data$allele.effects) |> as.data.frame()
  allele_effects_transposed$marker <- rownames(allele_effects_transposed)
  
  
  num_strains <- ncol(allele_effects_transposed) -1 
  chromosome_name <- paste0("chr", current_locus_info$chr)
  founder_strains <- rownames(current_scan_data$allele.effects)
  
  plot_data_list <- lapply(1:num_strains, function(strain){
    curr_strain <- colnames(allele_effects_transposed)[[strain]]
    temp_df <- allele_effects_transposed |> 
      dplyr::select(marker, founder_strains[strain]) %>% 
      filter(marker %in% markers$marker)
    
    temp_df <- left_join(temp_df, markers, by = "marker")
    
    temp_df <- temp_df |> 
      mutate(chrom = chromosome_name) |> 
      arrange(start)
    colnames(temp_df)[2] <- "score"
    
    temp_df$end <- c(temp_df$start[2:nrow(temp_df)] - 1L, 
                     temp_df$start[nrow(temp_df)] + 1)   
    temp_df <- temp_df |> 
      dplyr::select(chrom, start, end, score)
    return(temp_df)
  }
  )
  names(plot_data_list) <- founder_strains
  
  # Get the max abs haplotype effects to set y axis 
  max_allele_effect <- plot_data_list |> rbindlist() |> pull(score) |> abs() |> max() |> ceiling()
  
  strain_colors <- rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), length.out = num_strains)
  signalPlots <- c()
  
  signalPlots[[1]] <- plotgardener::plotSignal( #:length(plot_data_list)
    data = plot_data_list[[1]],
    params = params_genome,
    range = c(-max_allele_effect,max_allele_effect), 
    linecolor = strain_colors[1],
    fill = NA, 
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y + 2*PLOT_PARAMS$plot_height + 2*0.2,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), default.units = "inches",
    baseline = TRUE, 
    baseline.color = "grey")
  
  annoYaxis(plot = signalPlots[[1]], 
            at = c(-max_allele_effect, 0, max_allele_effect),
            axisLine     = TRUE,
            fontsize     = 8,
            main         = FALSE)
  lapply(2:8, function(strain_data){
    signalPlots[[strain_data]] <- plotgardener::plotSignal( #:length(plot_data_list)
      data = plot_data_list[[strain_data]],
      params = params_genome,
      range = c(-max_allele_effect,max_allele_effect), 
      linecolor = strain_colors[strain_data],
      fill = NA, 
      x              = PLOT_PARAMS$x,
      y              = PLOT_PARAMS$plot_y + 2*PLOT_PARAMS$plot_height + 2*0.2,
      width          = PLOT_PARAMS$plot_width,
      height         = PLOT_PARAMS$plot_height,
      just = c("center", "top"), default.units = "inches",
      baseline = TRUE, 
      baseline.color = "grey")
  }
  )

  plotText(label = "Founder Effects",  
           x       = 2 * PLOT_PARAMS$x + 0.1,
           y       = PLOT_PARAMS$plot_y + 2.5*PLOT_PARAMS$plot_height + 2*0.2,
           rot     = 270,
           fontsize= 8,
           just    = c("center","center"),
           default.units = "in")
  
  # --- Plot Genes ---
  nColors <- length(unique(overlaps$traitXdrug))
  ranges <- plotRanges(
    data = overlaps,
    params = params_genome,
    order = "random",
    fill = colorby("traitXdrug", 
                   palette = colorRampPalette(brewer.pal(nColors, name = "Set3"))),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y + 2.8*PLOT_PARAMS$plot_height + 3*0.2,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), default.units = "inches"
  )
  plotText(label = "Other Loci",  
           x       = 2 * PLOT_PARAMS$x + 0.1,
           y       = PLOT_PARAMS$plot_y + 3.5 *PLOT_PARAMS$plot_height + 3*0.2,
           rot     = 270,
           fontsize= 8,
           just    = "center",
           default.units = "in")

  # This requires `genes_in_locus` to be prepared
  genes_y_pos <- PLOT_PARAMS$plot_y + 4*PLOT_PARAMS$plot_height + 4*0.2 # Position below pyLMM
  founders <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/H1LtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")
  legendPlot <- plotLegend(
    legend = founders,
    fill = strain_colors,
    border = FALSE, 
    x = 8.7, y = PLOT_PARAMS$plot_y + 2.5*PLOT_PARAMS$plot_height + 2*0.2, width = 2, height = 1,
    fontsize = 7,
    just = c("left", "center"),
    orientation = "v",
    default.units = "inches"
  )
  
  legendPlot <- plotLegend(
    legend = unique(overlaps$traitXdrug),
    fill = brewer.pal(nColors, name = "Set3"),
    border = FALSE, 
    x = 8.7, y = PLOT_PARAMS$plot_y + 3.5*PLOT_PARAMS$plot_height + 3*0.2, width = 2, height = 1,
    fontsize = 7,
    just = c("left", "center"),
    orientation = "v",
    default.units = "inches"
  )
  legendPlot <- plotLegend(
    legend = c("Protein coding, human orthologue, and >5 cpm in NRVMs"),
    fill = "darkred",
    border = FALSE, 
    x = 0.25, y = genes_y_pos-0.1, width = 5, height = 0.15,
    fontsize = 7,
    just = c("left", "top"),
    orientation = "v",
    default.units = "inches"
  )
  
  gene_plot <- plotGenes(
    params = params_genome,
    x = 4.25, y = genes_y_pos, width = 8, height = 1, # Adjust height as needed
    just = c("center", "top"), default.units = "inches",
    geneOrder = top_genes_in_locus,
    fontsize = 6,
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
  trait_name <- current_locus_info$trait
  drug_condition <- current_locus_info$drug
  pheno_data = all_pheno_data
  output_dir = current_output_dir

  # Filter for the specific drug condition
  current_pheno_data <- pheno_data[Drug == drug_condition, ]
  
  plot_dir <- file.path(output_dir, "phenoDistributions")
  plot_file_name <- file.path(plot_dir, paste0("pheno_dist_", trait_name, "_", drug_condition, ".png"))
  
  if(!dir.exists(plot_dir)){
    dir.create(plot_dir, recursive = T)
  }
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




### Group Loci by chromsome 
sig_regions <- sig_regions |>
  mutate(start = pmin(upper_pos_lod_drop, lower_pos_lod_drop),
         end   = pmax(upper_pos_lod_drop, lower_pos_lod_drop))

gr <- GRanges(seqnames = paste0("chr", sig_regions$chr),        # add “chr” prefix so all seqlevels look alike
              ranges   = IRanges(sig_regions$start, sig_regions$end),
              mcols    = sig_regions)               # keeps all the original columns
hits <- GenomicRanges::findOverlaps(gr)                 # one‑based; 1, 2, 3, …

edges <- as.data.frame(hits)[ , 1:2]               # queryHits | subjectHits
edges <- subset(edges, queryHits != subjectHits)   # drop the trivial self‑edges
edges <- unique(edges)                             # remove duplicated pairs

g <- graph_from_data_frame(edges,
                           directed = FALSE,
                           vertices = data.frame(name = seq_along(gr)))

## 3. Connected components = overlapping loci --------------------------------
comp <- components(g)$membership                   # integer vector, length 42
sig_regions$locus_cluster <- comp

# --- 5. Main Loop to Generate Packets --- ####

# Define which loci to process.
group_to_process <- unique(sig_regions$locus_cluster)#[[2]] 

# Base output directory for all packets
base_output_dir <- "results/qtl_packets" 
if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

for (i in seq_along(group_to_process)) {
  i <- 2
  loci <- sig_regions |> filter(locus_cluster %in%   group_to_process[i])
  # Create a specific directory for this locus packet
  locus_packet_name <- paste0("locus_chr", loci$chr[[1]],
                              "_", round(min(loci$upper_pos_lod_drop), 0), "-",
                                   round(max(loci$lower_pos_lod_drop), 0), "Mb")
  current_output_dir <- file.path(base_output_dir, locus_packet_name)
  if (!dir.exists(current_output_dir)) dir.create(current_output_dir, recursive = TRUE)
  
  message(paste0("\n--- Processing Locus: ", locus_packet_name, " ---"))
  
  # --- A. Get genes and founder SNPs in the current locus region ---
  # Define locus boundaries (e.g., from LOD drop or fixed window around peak)
  # Using LOD drop from sig_regions, converting Mb to bp
  locus_start_bp <- min(loci$upper_pos_lod_drop) * 1e6
  locus_end_bp   <- max(loci$lower_pos_lod_drop)* 1e6
  
  genes_in_current_locus <- get_genes_in_region(loci$chr[[1]],
                                                locus_start_bp,
                                                locus_end_bp)
  
  founder_snps_current_locus <- get_founder_snps_in_region(loci$chr[[1]],
                                                           locus_start_bp,
                                                           locus_end_bp)
  
  top_genes_in_locus <- data.frame(gene = merged_gene_info |> 
    filter(mouse_gene_symbol %in% genes_in_current_locus$mouse_gene_symbol) |> 
    pull(mouse_gene_symbol),
  color = "#e34a33")
  
  # Make the overlaps dataframe to plot ranges of loci within loci cluster
  overlaps <- loci |> 
    mutate(chrom = paste0("chr", chr),
           start = upper_pos_lod_drop*10^6,
           end   = lower_pos_lod_drop*10^6, 
           strand= "-",
           traitXdrug = paste0(trait, ": ", drug)) |>  
    dplyr::select(chrom, start, end, strand, trait, drug, traitXdrug)
  
  # --- B. Generate Locus Zoom Plot ---
  # Constants
  COLORS <- c(Sig = "#1f78b4", nonSig = "#a6cee3")
  PLOT_DIMS <- list(page_width = 9, page_height = 6.5, res = 300)
  PLOT_PARAMS <- list(x = 4.25, plot_width = 8, plot_height = 1, plot_y = 0.5)
  GENE_DIMS <- list(height = 2, y_offset = 0.5, label_offset = 0.1)
  MIN_YLIM <- 5
  
  # --- B. Generate Locus Zoom Plot ---
  for(j in 1:nrow(loci)){
    print(j )
    current_locus_info <- loci[j,]
    locus_plot_file <- generate_locus_zoom_plot(current_locus_info,
                                                current_output_dir,
                                                genes_in_current_locus, 
                                                top_genes_in_locus, # Pass this for potential use in plotGenes
                                                founder_snps_current_locus) # Pass for founder variant track
  }
  
  # --- C. Generate Phenotype Distribution Plot ---
  for(j in 1:nrow(loci)){
    current_locus_info <- loci[j,]
  pheno_dist_plot_file <- generate_phenotype_dist_plot(current_locus_info$trait,
                                                       current_locus_info$drug,
                                                       pheno_data = all_pheno_data,
                                                       output_dir = current_output_dir)
  }
  # --- D. to generate Excel file
  # --- D. Generate the combined gene info Excel file for all loci in the cluster ---
  # Pass the entire loci dataframe instead of a single locus
  gene_info_excel_file <- generate_gene_info_excel(
    genes_in_current_locus,
    loci_info_df = loci,  # Pass the entire loci dataframe 
    current_output_dir,
    full_mouse_pheno_data = mouse_pheno,
    full_human_assoc_data = associations,
    orthology_data = ortho_mouse2h
  )
  # --- E. Generate Founder Mutation Table for SNPs in Locus ---
  founder_snp_table_file <- generate_founder_mutation_table(founder_snps_current_locus,
                                                            current_locus_info,
                                                            current_output_dir)
  

  # Define heart-related terms for human disease associations
  human_heart_terms <- c(
    # Basic cardiac terms
    "heart", "cardi", "coronary", "atri", "ventric", "myocard", 
    # Heart function issues
    "hypertroph", "fibril", "valv", "rhythm", "arrhythm", "tachycard", "bradycard",
    # Circulation and vascular
    "blood pressure", "pulse", "aort", "angina", "stroke", "hypertens", "thromb", "pulse pressure",
    # Ischemia and infarction
    "infarct", "ischem"
  )
  
  # Define heart-related MGI phenotype terms with more comprehensive coverage
  mgi_heart_terms <- c(
    # Heart structure terms
    "heart weight", "heart atrium", "heart ventricle", "heart morphology", "myocardi", 
    "cardiac", "cardio", "atrioventricular", "pericardial", "ventricular septal defect",
    "atrial", "thin myocardium", "common atrioventricular valve",
    
    # Vasculature terms
    "coronary vessel", "blood vessel", "vascular", "vasculature", "aort",
    "lymphatic vessel", "yolk sac vascular", "arch artery", "pharyngeal arch artery",
    
    # Circulation terms
    "heart failure", "circulat", "thromb", "hemorrhage", "congest", "blood circulation",
    
    # Heart function terms
    "response of heart", "cardiac muscle", "myocardium layer", "myocardial fiber",
    
    # General pattern matches for common phenotype descriptions
    "dilated heart", "increased.+heart", "decreased.+heart", "abnormal.+heart", 
    "common.+valve", "interrupted aortic arch"
  )
  
  # Create patterns - using word boundaries for more precise matches
  human_heart_pattern <- paste0("\\b(", paste(human_heart_terms, collapse="|"), ")\\b|", "LV\\.")
  mgi_heart_pattern <- paste0("\\b(", paste(mgi_heart_terms, collapse="|"), ")\\b")
  
  # Filter for human heart/cardiology-related genes - carefully considering what constitutes a match
  human_genes <- gene_info_excel_file %>%
    filter(grepl(human_heart_pattern, `Human Disease Summary`, ignore.case = TRUE))
  
  # Filter for mouse heart/cardiology-related genes
  mouse_genes <- gene_info_excel_file %>%
    filter(grepl(mgi_heart_pattern, `Mouse Phenotype Summary (MGI)`, ignore.case = TRUE))
  
  # Get subsets of genes present in one or both
  mouse_only_genes <- mouse_genes %>% 
    filter(!`Mouse Gene Symbol` %in% human_genes$`Mouse Gene Symbol`) %>%
    arrange(`Mouse Gene Symbol`)
  
  human_only_genes <- human_genes %>% 
    filter(!`Mouse Gene Symbol` %in% mouse_genes$`Mouse Gene Symbol`) %>%
    arrange(`Mouse Gene Symbol`)
  
  human_mouse_genes <- human_genes %>% 
    filter(`Mouse Gene Symbol` %in% mouse_genes$`Mouse Gene Symbol`) %>%
    arrange(`Mouse Gene Symbol`)
  
  # Function to extract cardiac traits from Human Disease Summary with careful parsing
  extract_human_cardiac_traits <- function(disease_summary) {
    if(is.na(disease_summary) || disease_summary == "") return("No cardiac traits identified")
    
    # Split the disease summary by commas and trim whitespace
    all_traits <- trimws(unlist(strsplit(disease_summary, ",")))
    
    # Filter for cardiac traits with more precise pattern matching
    cardiac_traits <- all_traits[grepl(human_heart_pattern, all_traits, ignore.case = TRUE)]
    
    # Return unique traits to avoid repetition
    return(paste(unique(cardiac_traits), collapse=", "))
  }
  
  # Function to extract cardiac phenotypes from Mouse Phenotype Summary with precise matching
  extract_mouse_cardiac_traits <- function(phenotype_summary) {
    if(is.na(phenotype_summary) || phenotype_summary == "") return("No cardiac phenotypes identified")
    
    # Split the phenotype summary by commas and trim whitespace
    all_phenotypes <- trimws(unlist(strsplit(phenotype_summary, ", ")))
    
    # Filter for cardiac phenotypes
    cardiac_phenotypes <- all_phenotypes[grepl(mgi_heart_pattern, all_phenotypes, ignore.case = TRUE)]
    
    # Return unique phenotypes to avoid repetition
    return(paste(unique(cardiac_phenotypes), collapse=", "))
  }
  
  # Format human-only genes and traits with careful spacing and formatting
  format_human_only <- function(genes_df) {
    if(nrow(genes_df) == 0) return("None identified")
    
    result <- character(nrow(genes_df))
    for(i in 1:nrow(genes_df)) {
      gene <- genes_df$`Mouse Gene Symbol`[i]
      traits <- extract_human_cardiac_traits(genes_df$`Human Disease Summary`[i])
      result[i] <- paste0(gene, ": ", traits)
    }
    
    return(paste(result, collapse="\n        "))
  }
  
  # Format mouse-only genes and phenotypes with consistent formatting
  format_mouse_only <- function(genes_df) {
    if(nrow(genes_df) == 0) return("None identified")
    
    result <- character(nrow(genes_df))
    for(i in 1:nrow(genes_df)) {
      gene <- genes_df$`Mouse Gene Symbol`[i]
      phenotypes <- extract_mouse_cardiac_traits(genes_df$`Mouse Phenotype Summary (MGI)`[i])
      result[i] <- paste0(gene, ": ", phenotypes)
    }
    
    return(paste(result, collapse="\n        "))
  }
  
  # Format genes with both human and mouse cardiac associations with hierarchical organization
  format_both <- function(genes_df) {
    if(nrow(genes_df) == 0) return("None identified")
    
    result <- character(nrow(genes_df))
    for(i in 1:nrow(genes_df)) {
      gene <- genes_df$`Mouse Gene Symbol`[i]
      human_traits <- extract_human_cardiac_traits(genes_df$`Human Disease Summary`[i])
      mouse_phenotypes <- extract_mouse_cardiac_traits(genes_df$`Mouse Phenotype Summary (MGI)`[i])
      result[i] <- paste0(gene, ":\n            Human: ", human_traits, "\n            Mouse: ", mouse_phenotypes)
    }
    
    return(paste(result, collapse="\n\n        "))
  }
  
  # Add summary statistics for better context
  human_only_count <- nrow(human_only_genes)
  mouse_only_count <- nrow(mouse_only_genes)
  both_count <- nrow(human_mouse_genes)
  total_count <- human_only_count + mouse_only_count + both_count
  
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
    "\n- All 'Associated Traits (Drug)' in the original data are cardiac-related QTL mappings",
    "\n- LV-related measurements refer to left ventricular parameters",
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
        Expression_Level = human_only_genes$`Avg NRVM CPM (Ctrl)`,
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
        Expression_Level = mouse_only_genes$`Avg NRVM CPM (Ctrl)`,
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
        Expression_Level = human_mouse_genes$`Avg NRVM CPM (Ctrl)`,
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
}

message("\nAll selected loci processed. Packets generated in: ", base_output_dir)

