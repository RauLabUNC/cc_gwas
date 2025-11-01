#' PackRat: Locus Visualization Functions
#'
#' LocusZoom-style plotting functions using plotgardener. These functions are
#' designed for mouse QTL data but can be adapted for other organisms.
#'
#' @author Brian Gural, Anh Luu, Todd Kimball, Christoph Rau
#' @references Gural et al. (in preparation)

library(plotgardener)
library(tidyverse)
library(GenomicRanges)
library(RColorBrewer)

# =============================================================================
# Main Locus Plot Function
# =============================================================================

#' Generate LocusZoom-style plot for a QTL locus
#'
#' Creates a multi-panel plot showing QTL scan statistics, founder allele effects,
#' overlapping loci, and gene tracks.
#'
#' @param locus_info List or data frame row with locus information
#'   Required fields: chr, peak_pos, upper_pos_lod_drop, lower_pos_lod_drop, trait
#' @param scan_data List with QTL scan results (LOD scores, positions, allele effects)
#' @param threshold_data List with significance thresholds
#' @param overlapping_loci Optional GRanges or data frame with other loci in region
#' @param genes_to_highlight Optional character vector of gene symbols to highlight
#' @param output_file Character, path to output PDF file
#' @param genome_assembly plotgardener assembly object (default mm39)
#' @param plot_width Numeric, plot width in inches (default 10.5)
#' @param plot_height Numeric, plot height in inches (default 5.5)
#' @param padding_bp Numeric, bp to add on either side of locus (default 500kb)
#'
#' @return Invisibly returns the output file path
#'
#' @details
#' This function creates a publication-quality LocusZoom plot with:
#' - Top panel: LOD scores across the region with significance threshold
#' - Middle panel: Founder strain allele effects
#' - Lower panels: Overlapping QTL regions and gene annotations
#'
#' Requires plotgardener, TxDb, and OrgDb packages for the organism.
#'
#' @export
plot_locus_zoom <- function(locus_info,
                           scan_data,
                           threshold_data,
                           overlapping_loci = NULL,
                           genes_to_highlight = NULL,
                           output_file,
                           genome_assembly = NULL,
                           plot_width = 10.5,
                           plot_height = 5.5,
                           padding_bp = 500000) {

  # Set default assembly if not provided
  if (is.null(genome_assembly)) {
    genome_assembly <- assembly(
      Genome = "mm39_GRCm39",
      TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene",
      OrgDb = "org.Mm.eg.db"
    )
  }

  # Extract locus parameters
  chr <- locus_info$chr
  peak_pos_mb <- locus_info$peak_pos
  upper_mb <- locus_info$upper_pos_lod_drop
  lower_mb <- locus_info$lower_pos_lod_drop
  trait <- locus_info$trait

  # Convert to bp and add padding
  plot_start_bp <- max(0, floor(upper_mb * 1e6) - padding_bp)
  plot_end_bp <- ceiling(lower_mb * 1e6) + padding_bp
  bounds_bp <- c(upper_mb * 1e6, lower_mb * 1e6)

  # Prepare plot file
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  message("Generating locus zoom plot: ", output_file)

  # --- Start plotting ---
  pdf(output_file, width = plot_width, height = plot_height)
  pageCreate(width = plot_width, height = plot_height,
            default.units = "inches", showGuides = FALSE)

  # Set up genome parameters
  params_genome <- pgParams(
    assembly = genome_assembly,
    chrom = paste0("chr", chr),
    chromstart = plot_start_bp,
    chromend = plot_end_bp
  )

  # --- Plot LOD scores ---
  lod_plot <- plot_lod_panel(
    scan_data = scan_data,
    threshold = threshold_data,
    params_genome = params_genome,
    chr = chr,
    bounds_bp = bounds_bp,
    y_position = 0.5
  )

  # --- Plot founder allele effects ---
  effects_plot <- plot_founder_effects(
    scan_data = scan_data,
    params_genome = params_genome,
    chr = chr,
    plot_start_bp = plot_start_bp,
    plot_end_bp = plot_end_bp,
    y_position = 0.5 + 1.2
  )

  # --- Plot overlapping loci (if provided) ---
  if (!is.null(overlapping_loci)) {
    overlaps_plot <- plot_overlapping_loci(
      overlapping_loci = overlapping_loci,
      params_genome = params_genome,
      y_position = 0.5 + 2.4
    )
  }

  # --- Plot genes ---
  genes_plot <- plot_gene_track(
    params_genome = params_genome,
    genes_to_highlight = genes_to_highlight,
    y_position = 0.5 + 3.6
  )

  # --- Add genome axis ---
  annoGenomeLabel(
    plot = genes_plot,
    params = params_genome,
    x = 4.25, y = 0.5 + 4.6,
    scale = "Mb", fontsize = 10,
    just = c("center", "top"),
    default.units = "inches"
  )

  # --- Add title ---
  plotText(
    label = paste("QTL:", trait, "- Chr", chr,
                 "Peak:", round(peak_pos_mb, 2), "Mb"),
    x = plot_width / 2, y = 0.1,
    just = c("center", "top"),
    fontface = "bold", fontsize = 12,
    default.units = "inches"
  )

  dev.off()

  invisible(output_file)
}


# =============================================================================
# Panel-Specific Plotting Functions
# =============================================================================

#' Plot LOD score panel
#' @keywords internal
plot_lod_panel <- function(scan_data, threshold, params_genome, chr, bounds_bp,
                          y_position = 0.5, plot_width = 8, plot_height = 1) {

  # Prepare LOD data for plotting
  lod_df <- tibble(
    marker = names(scan_data$LOD),
    chr = as.character(scan_data$chr),
    pos = scan_data$pos$Mb * 1e6,
    lod = scan_data$LOD
  ) %>%
    filter(chr == !!chr, !is.na(pos), !is.na(lod)) %>%
    transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod))

  # Determine y-axis limits
  ylim <- c(0, max(c(-log10(lod_df$p), threshold, 5), na.rm = TRUE) + 1)

  # Plot Manhattan
  lod_plot <- plotManhattan(
    data = lod_df,
    params = params_genome,
    range = ylim,
    trans = "-log10",
    sigVal = 10^(-threshold),
    x = 4.25, y = y_position,
    width = plot_width, height = plot_height,
    just = c("center", "top"),
    xfield = "pos", yfield = "p",
    fill = "#a6cee3", sigCol = "#1f78b4",
    sigLine = TRUE, baseline = TRUE,
    default.units = "inches"
  )

  # Add y-axis
  annoYaxis(
    plot = lod_plot,
    at = pretty(ylim),
    axisLine = TRUE,
    fontsize = 8,
    main = FALSE
  )

  # Add y-axis label
  plotText(
    label = "LOD",
    x = 0.3, y = y_position + plot_height / 2,
    rot = 270, fontsize = 8,
    just = "center",
    default.units = "inches"
  )

  # Highlight significant region
  annoHighlight(
    plot = lod_plot,
    chrom = paste0("chr", chr),
    chromstart = floor(min(bounds_bp)),
    chromend = ceiling(max(bounds_bp)),
    fill = "#fb9a99",
    y = y_position,
    height = plot_height,
    just = c("left", "top"),
    default.units = "inches",
    alpha = 0.2,
    params = params_genome
  )

  return(lod_plot)
}


#' Plot founder allele effects panel
#' @keywords internal
plot_founder_effects <- function(scan_data, params_genome, chr,
                                 plot_start_bp, plot_end_bp,
                                 y_position = 1.7, plot_width = 8, plot_height = 1) {

  # Get markers in range
  in_range <- scan_data$chr == chr &
              scan_data$pos$Mb * 1e6 >= plot_start_bp &
              scan_data$pos$Mb * 1e6 <= plot_end_bp

  markers <- data.frame(
    marker = scan_data$loci[in_range],
    start = scan_data$pos$Mb[in_range] * 1e6
  ) %>% arrange(start)

  # Prepare founder effects data
  founder_strains <- rownames(scan_data$allele.effects)
  n_strains <- length(founder_strains)

  plot_data_list <- lapply(1:n_strains, function(i) {
    strain_effects <- data.frame(
      marker = colnames(scan_data$allele.effects),
      score = scan_data$allele.effects[i, ]
    ) %>%
      filter(marker %in% markers$marker) %>%
      left_join(markers, by = "marker") %>%
      mutate(chrom = paste0("chr", chr)) %>%
      arrange(start)

    # Create ranges
    strain_effects$end <- c(strain_effects$start[2:nrow(strain_effects)] - 1,
                           strain_effects$start[nrow(strain_effects)] + 1)

    strain_effects %>% select(chrom, start, end, score)
  })
  names(plot_data_list) <- founder_strains

  # Get max effect for y-axis
  max_effect <- plot_data_list %>%
    bind_rows() %>%
    pull(score) %>%
    abs() %>%
    max(na.rm = TRUE) %>%
    {round(. * 1.1, 2)}

  # Founder colors
  strain_colors <- rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                        "#66A61E", "#E6AB02", "#A6761D", "#666666"),
                      length.out = n_strains)

  # Plot first strain
  effects_plot <- plotSignal(
    data = plot_data_list[[1]],
    params = params_genome,
    range = c(-max_effect, max_effect),
    linecolor = strain_colors[1],
    fill = NA,
    x = 4.25, y = y_position,
    width = plot_width, height = plot_height,
    just = c("center", "top"),
    default.units = "inches",
    baseline = TRUE,
    baseline.color = "grey"
  )

  # Overlay remaining strains
  for (i in 2:n_strains) {
    plotSignal(
      data = plot_data_list[[i]],
      params = params_genome,
      range = c(-max_effect, max_effect),
      linecolor = strain_colors[i],
      fill = NA,
      x = 4.25, y = y_position,
      width = plot_width, height = plot_height,
      just = c("center", "top"),
      default.units = "inches",
      baseline = TRUE,
      baseline.color = "grey"
    )
  }

  # Add y-axis
  annoYaxis(
    plot = effects_plot,
    at = c(-max_effect, 0, max_effect),
    axisLine = TRUE,
    fontsize = 8,
    main = FALSE
  )

  # Add y-axis label
  plotText(
    label = "Founder Effects",
    x = 0.3, y = y_position + plot_height / 2,
    rot = 270, fontsize = 8,
    just = "center",
    default.units = "inches"
  )

  # Add legend
  plotLegend(
    legend = founder_strains,
    fill = strain_colors,
    border = FALSE,
    x = 8.7, y = y_position,
    width = 2, height = 1,
    fontsize = 7,
    just = c("left", "center"),
    orientation = "v",
    default.units = "inches"
  )

  return(effects_plot)
}


#' Plot overlapping QTL loci
#' @keywords internal
plot_overlapping_loci <- function(overlapping_loci, params_genome,
                                  y_position = 2.9, plot_width = 8, plot_height = 1) {

  # Ensure overlapping_loci is in proper format
  if (!is(overlapping_loci, "GRanges")) {
    overlapping_loci <- GRanges(
      seqnames = paste0("chr", overlapping_loci$chr),
      ranges = IRanges(overlapping_loci$start, overlapping_loci$end),
      trait = overlapping_loci$trait
    )
  }

  n_loci <- length(unique(overlapping_loci$trait))
  loci_colors <- colorRampPalette(brewer.pal(min(n_loci, 12), "Set3"))(n_loci)

  ranges_plot <- plotRanges(
    data = overlapping_loci,
    params = params_genome,
    order = "random",
    fill = colorby("trait", palette = loci_colors),
    x = 4.25, y = y_position,
    width = plot_width, height = plot_height,
    just = c("center", "top"),
    default.units = "inches"
  )

  # Add y-axis label
  plotText(
    label = "Other QTL",
    x = 0.3, y = y_position + plot_height / 2,
    rot = 270, fontsize = 8,
    just = "center",
    default.units = "inches"
  )

  return(ranges_plot)
}


#' Plot gene track with highlighting
#' @keywords internal
plot_gene_track <- function(params_genome, genes_to_highlight = NULL,
                           y_position = 4.1, plot_width = 8, plot_height = 1) {

  # Prepare highlight data frame
  if (!is.null(genes_to_highlight)) {
    gene_highlight_df <- data.frame(
      gene = genes_to_highlight,
      color = "#e34a33"
    )
  } else {
    gene_highlight_df <- NULL
  }

  gene_plot <- plotGenes(
    params = params_genome,
    x = 4.25, y = y_position,
    width = plot_width, height = plot_height,
    just = c("center", "top"),
    default.units = "inches",
    geneOrder = if (!is.null(genes_to_highlight)) genes_to_highlight else NULL,
    fontsize = 6,
    geneHighlights = if (!is.null(genes_to_highlight)) gene_highlight_df else NULL,
    geneBackground = "#fdbb84"
  )

  return(gene_plot)
}
