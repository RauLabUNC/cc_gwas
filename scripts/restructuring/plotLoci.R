# QTL Plotting Script
# Compares results of BoxCox and Z score norm for trait data

# Load Libraries
library(tidyverse)
library(plotgardener)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)

# Constants
COLORS <- list(
  boxcox = c(Sig = "#1f78b4", nonSig = "#a6cee3"),
  zscore = c(Sig = "#ff7f00", nonSig = "#fdbf6f")
)
PLOT_DIMS <- list(page_width = 7, page_height = 7, res = 300)
PLOT_PARAMS <- list(x = 3.25, plot_width = 6, plot_height = 2.2, plot_y = 0.5)
GENE_DIMS <- list(height = 1.5, y_offset = 0.25, label_offset = 0.1)
MIN_YLIM <- 5

# Helper: Compute padding
compute_padding <- function(start, end, min_pad = 5e5, pad_frac = 0.1) {
  width <- max(1, end - start)
  pad <- max(min_pad, width * pad_frac)
  c(max(0, floor(start - pad)), ceiling(end + pad))
}

# Safe loaders using purrr
get_scan_safe <- purrr::safely(
  function(norm, trait, drug) {
    scan_raw <- scan_data[[norm]][[trait]][[drug]]
    tibble(
      marker = names(scan_raw$LOD),
      chr    = as.character(scan_raw$chr),
      pos    = scan_raw$pos$Mb * 1e6,
      lod    = scan_raw$LOD
    ) %>%
      filter(!is.na(chr) & !is.na(pos) & !is.na(lod))
  }
)
get_thresh_safe <- purrr::safely(
  function(norm, trait, drug) {
    threshold_data[[norm]][[trait]][[drug]]
  }
)

# Load Data
sig_regions <- read_csv("results/sig_regions/all_significant_regions_summary.csv") %>%
  mutate(
    across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
    chr = as.character(chr)
  ) %>%
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, norm, drug)
scan_data      <- readRDS("results/sig_regions/scan_data.rds")
threshold_data <- readRDS("results/sig_regions/threshold_data.rds")

# Genome assembly
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

# Output directory
output_dir <- "results/sig_loci_plots_comparison_v3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Plotting function
add_lod_curve <- function(params_genome, scan_df, norm, threshold, type = "orig") {
  df <- scan_df %>%
    transmute(
      chrom = paste0("chr", chr),
      pos, p = 10^(-lod)
    )
  y_origin <- PLOT_PARAMS$plot_y + ifelse(type == "alt", PLOT_PARAMS$plot_height + 0.2, 0)
  plotManhattan(
    data = df,
    params   = params_genome,
    trans    = "-log10",
    sigVal   = 10^(-threshold),
    range    = plot_ylim,
    x        = PLOT_PARAMS$x,
    y        = y_origin,
    width    = PLOT_PARAMS$plot_width,
    height   = PLOT_PARAMS$plot_height,
    just     = c("center", "top"),
    xfield   = "pos",
    yfield   = "p",
    fill     = COLORS[[norm]]["nonSig"],
    sigCol   = COLORS[[norm]]["Sig"],
    sigLine  = TRUE,
    baseline = TRUE,
    yscale_reverse = FALSE
  )
}

# Main loop over significant regions
sig_regions %>%
  purrr::pwalk(function(chr, norm, trait, drug, peak_pos, lower_pos_lod_drop, upper_pos_lod_drop, ...) {
    # Prepare chromosome string and bounds
    chr_str <- if (!startsWith(chr, "chr")) paste0("chr", chr) else chr
    bounds_bp <- c(lower_pos_lod_drop * 1e6, upper_pos_lod_drop * 1e6)
    padded    <- compute_padding(min(bounds_bp), max(bounds_bp))
    chromstart <- padded[1]
    chromend   <- padded[2]
    
    # Genomic plotting parameters
    params_genome <- pgParams(
      assembly   = mm39,
      chrom      = chr_str,
      chromstart = chromstart,
      chromend   = chromend
    )
    
    # Alternative conditions
    alt_drug <- if (drug == "Ctrl") "Iso" else "Ctrl"
    alt_norm <- if (norm  == "boxcox") "zscore" else "boxcox"
    
    # Retrieve scan & threshold data
    orig_scan       <- get_scan_safe(norm, trait, drug)$result
    scan_alt_drug   <- get_scan_safe(norm, trait, alt_drug)$result
    scan_alt_norm   <- get_scan_safe(alt_norm, trait, drug)$result
    scan_alt_both   <- get_scan_safe(alt_norm, trait, alt_drug)$result
    
    thresh_orig       <- get_thresh_safe(norm, trait, drug)$result
    thresh_alt_drug   <- get_thresh_safe(norm, trait, alt_drug)$result
    thresh_alt_norm   <- get_thresh_safe(alt_norm, trait, drug)$result
    thresh_alt_both   <- get_thresh_safe(alt_norm, trait, alt_drug)$result
    
    # Determine y-axis limits
    lods   <- c(orig_scan$lod, scan_alt_drug$lod, scan_alt_norm$lod, scan_alt_both$lod)
    thresh <- c(thresh_orig, thresh_alt_drug, thresh_alt_norm, thresh_alt_both)
    plot_ylim <<- c(
      0,
      max(MIN_YLIM, ceiling(max(c(lods, thresh), na.rm = TRUE)))
    )
    
    # Output filename
    outfile <- file.path(
      output_dir,
      paste0("locus_", norm, "_", trait, "_", drug,
             "_chr", chr, "_peak", round(peak_pos, 2), "Mb.png")
    )
    
    # Start plotting
    png(outfile, width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height,
        units = "in", res = PLOT_DIMS$res)
    pageCreate(
      width = PLOT_DIMS$page_width,
      height = PLOT_DIMS$page_height,
      default.units = "inches",
      showGuides    = FALSE
    )
    
    # Plot curves
    main_plot <- add_lod_curve(params_genome, orig_scan, norm, thresh_orig, type = "orig")
    annoYaxis(plot    = main_plot,
              at      = pretty(plot_ylim),
              axisLine= TRUE,
              fontsize= 8,
              main    = FALSE)
    
    add_lod_curve(params_genome, scan_alt_norm, norm = alt_norm, threshold = thresh_alt_norm)
    sub_plot <- add_lod_curve(params_genome, scan_alt_both, norm = alt_norm, threshold = thresh_alt_both, type = "alt")
    add_lod_curve(params_genome, scan_alt_drug, norm = norm, threshold = thresh_alt_drug, type = "alt")
    annoYaxis(plot    = sub_plot,
              at      = pretty(plot_ylim),
              axisLine= TRUE,
              fontsize= 8,
              main    = FALSE)
    
    # Annotations & highlights
    plotText(label   = paste0("LOD — ", drug),
             x       = 2 * PLOT_PARAMS$x + 0.1,
             y       = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height / 2,
             rot     = 270,
             fontsize= 10,
             just    = "center",
             default.units = "in")
    plotText(label   = paste0("LOD — ", alt_drug),
             x       = 2 * PLOT_PARAMS$x + 0.1,
             y       = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height / 2 + PLOT_PARAMS$plot_height + 0.2,
             rot     = 270,
             fontsize= 10,
             just    = "center",
             default.units = "in")
    
    annoHighlight(
      plot       = main_plot,
      chrom      = chr_str,
      chromstart = floor(min(bounds_bp)),
      chromend   = ceiling(max(bounds_bp)),
      fill       = "#fb9a99",
      y          = PLOT_PARAMS$plot_y,
      height     = PLOT_PARAMS$plot_height,
      just       = c("left", "top"),
      default.units = "in",
      alpha      = 0.2,
      params     = params_genome
    )
    
    genes_plot <- plotGenes(
      params        = params_genome,
      x             = PLOT_PARAMS$x,
      y             = PLOT_PARAMS$plot_y + 2*PLOT_PARAMS$plot_height + GENE_DIMS$y_offset,
      width         = PLOT_PARAMS$plot_width,
      height        = GENE_DIMS$height,
      stroke        = 1,
      fontsize      = 10,
      just          = c("center", "top"),
      default.units = "in"
    )
    
    annoGenomeLabel(
      plot        = genes_plot,
      params      = params_genome,
      x           = PLOT_PARAMS$x,
      y           = PLOT_PARAMS$plot_y + 2*PLOT_PARAMS$plot_height + GENE_DIMS$height + GENE_DIMS$y_offset,
      scale       = "Mb",
      fontsize    = 8,
      just        = c("center", "top"),
      default.units = "in"
    )
    
    plotText(
      label   = paste("QTL:", trait, "(", norm, "/", drug, ") -", chr_str,
                      "Peak:", round(peak_pos, 2), "Mb"),
      x       = PLOT_PARAMS$x,
      y       = 0.05,
      just    = c("center", "top"),
      fontface= "bold",
      fontsize= 12,
      default.units = "in"
    )
    
    plotLegend(
      legend    = c("Box-Cox", "Z-score"),
      fill      = c(COLORS$boxcox["Sig"], COLORS$zscore["Sig"]),
      border    = FALSE,
      x         = PLOT_PARAMS$x,
      y         = -0.1,
      width     = 1,
      height    = 0.75,
      just      = c("center", "top"),
      fontsize  = 10,
      default.units = "in",
      orientation= "h"
    )
    
    dev.off()
  })