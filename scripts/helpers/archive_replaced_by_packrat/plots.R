

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
