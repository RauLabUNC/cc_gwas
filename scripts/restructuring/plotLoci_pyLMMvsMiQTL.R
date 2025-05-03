# Refactored QTL Plotting Script with pyResults Sub-Plot

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
PLOT_DIMS <- list(page_width = 7, page_height = 6.5, res = 300)
PLOT_PARAMS <- list(x = 3.25, plot_width = 6, plot_height = 1.5, plot_y = 0.5)
GENE_DIMS <- list(height = 2, y_offset = 0.5, label_offset = 0.1)
MIN_YLIM <- 5

# Helper: Compute padding around region
compute_padding <- function(start, end, min_pad = 5e5, pad_frac = 0.1) {
  width <- max(1, end - start)
  pad   <- max(min_pad, width * pad_frac)
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
sig_regions    <- read_csv("results/sig_regions/all_significant_regions_summary.csv") %>%
  mutate(
    across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
    chr = as.character(chr)
  ) %>%
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, norm, drug)
scan_data      <- readRDS("results/sig_regions/scan_data.rds")
threshold_data <- readRDS("results/sig_regions/threshold_data.rds")


# pyLMM results
pyResults <- c()

pyResults[["Ctrl"]] <- read.csv("data/processed/christophGWAS/Ctrl_pvals.csv")
pyResults[["Iso"]] <- read.csv("data/processed/christophGWAS/Iso_pvals.csv")


# Genome assembly
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

# Ensure pyResults is in environment
# pyResults is a list: pyResults$Ctrl, pyResults$Iso data.frames

# Output directory
output_dir <- "results/miQTL_pyLMM"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Plotting function for original LOD curves
add_lod_curve <- function(params_genome, scan_df, norm, threshold, type = "orig") {
  df <- scan_df %>%
    transmute(
      chrom = paste0("chr", chr),
      pos, p = 10^(-lod)
    )
  y_origin <- PLOT_PARAMS$plot_y + ifelse(type == "alt", PLOT_PARAMS$plot_height + 0.2, 0)
  plotManhattan(
    data           = df,
    params         = params_genome,
    trans          = "-log10",
    sigVal         = 10^(-threshold),
    range          = plot_ylim,
    x              = PLOT_PARAMS$x,
    y              = y_origin,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just           = c("center", "top"),
    xfield         = "pos",
    yfield         = "p",
    fill           = COLORS[[norm]]["nonSig"],
    sigCol         = COLORS[[norm]]["Sig"],
    sigLine        = TRUE,
    baseline       = TRUE,
    yscale_reverse = FALSE,
    default.units  = "inches"
  )
}
# Main loop over significant regions
sig_regions %>%
  purrr::pwalk(function(chr, norm, trait, drug, peak_pos, lower_pos_lod_drop, upper_pos_lod_drop, ...) {
    # Chromosome string and bounds
    chr_str   <- if (!startsWith(chr, "chr")) paste0("chr", chr) else chr
    bounds_bp <- c(lower_pos_lod_drop * 1e6, upper_pos_lod_drop * 1e6)
    padded    <- compute_padding(min(bounds_bp), max(bounds_bp))
    chromstart <- padded[1]
    chromend   <- padded[2]
    
    # Genomic params
    params_genome <- pgParams(
      assembly   = mm39,
      chrom      = chr_str,
      chromstart = chromstart,
      chromend   = chromend
    )
    
    # Alternative conditions
    alt_drug <- if (drug == "Ctrl") "Iso" else "Ctrl"
    alt_norm <- if (norm  == "boxcox") "zscore" else "boxcox"
    
    # Retrieve original scan & thresholds
    orig_scan     <- get_scan_safe(norm, trait, drug)$result
    thresh_orig   <- get_thresh_safe(norm, trait, drug)$result
    
    # Define y-axis limits for main plot
    lods_main <- orig_scan$lod
    plot_ylim <- c(
      0,
      max(MIN_YLIM, ceiling(max(lods_main, thresh_orig, na.rm = TRUE)+0.5))
    )
    
    # Prepare output file
    outfile <- file.path(
      output_dir,
      paste0("locus_", norm, "_", trait, "_", drug,
             "_chr", chr, "_peak", round(peak_pos, 2), "Mb.png")
    )
    
    # Begin plotting
    png(outfile,
        width       = PLOT_DIMS$page_width,
        height      = PLOT_DIMS$page_height,
        units       = "in",
        res         = PLOT_DIMS$res)
    pageCreate(
      width         = PLOT_DIMS$page_width,
      height        = PLOT_DIMS$page_height,
      default.units = "inches",
      showGuides    = FALSE
    )
    
    # Main LOD curve
    main_plot <- add_lod_curve(params_genome, orig_scan, norm, thresh_orig, type = "orig")
    annoYaxis(
      plot         = main_plot,
      at           = pretty(plot_ylim),
      axisLine     = TRUE,
      fontsize     = 8,
      main         = FALSE
    )
    # Annotations & highlights
    plotText(label   = paste0("LOD, miQTL — ", drug),
             x       = 2 * PLOT_PARAMS$x + 0.1,
             y       = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height / 2,
             rot     = 270,
             fontsize= 8,
             just    = "center",
             default.units = "in")
    
    # Sub-plot: use pyResults (no threshold)
    py_df <- pyResults[[drug]] %>%
      filter(Chr == chr) %>%
      transmute(
        chrom = paste0("chr", Chr),
        pos   = Pos_mm39,
        p = .data[[trait]]
      )
    # y-axis limits for sub-plot
    sub_ylim <- c(
      0,
      ceiling(max(-log10(py_df$p), na.rm = TRUE))
    )
    sub_plot <- plotManhattan(
      data           = py_df,
      params         = params_genome,
      trans          = "-log10",
      range          = sub_ylim,
      x              = PLOT_PARAMS$x,
      y              = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height + 0.2,
      width          = PLOT_PARAMS$plot_width,
      height         = PLOT_PARAMS$plot_height,
      just           = c("center", "top"),
      xfield         = "pos",
      yfield         = "value",
      fill           = "black",
      baseline       = TRUE,
      sigLine        = FALSE,
      default.units  = "inches"
    )
    annoYaxis(
      plot         = sub_plot,
      at           = pretty(sub_ylim),
      axisLine     = TRUE,
      fontsize     = 8,
      main         = FALSE
    )
    plotText(label   = paste0("-log10(p), pyLMM — ", drug),
             x       = 2 * PLOT_PARAMS$x + 0.1,
             y       = PLOT_PARAMS$plot_y + 1.5*PLOT_PARAMS$plot_height + 0.2,
             rot     = 270,
             fontsize= 8,
             just    = "center",
             default.units = "in")
    
    # Highlight significant region on main_plot
    annoHighlight(
      plot         = main_plot,
      chrom        = chr_str,
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
    
    # Genes and genome label
    genes_plot <- plotGenes(
      params        = params_genome,
      x             = PLOT_PARAMS$x,
      y             = PLOT_PARAMS$plot_y + 2 * PLOT_PARAMS$plot_height + GENE_DIMS$y_offset,
      width         = PLOT_PARAMS$plot_width,
      height        = GENE_DIMS$height,
      stroke        = 1,
      fontsize      = 10,
      just          = c("center", "top"),
      default.units = "inches"
    )
    annoGenomeLabel(
      plot        = genes_plot,
      params      = params_genome,
      x           = PLOT_PARAMS$x,
      y           =  PLOT_PARAMS$plot_y + 2 * PLOT_PARAMS$plot_height + GENE_DIMS$y_offset + GENE_DIMS$height,
      scale       = "Mb",
      fontsize    = 8,
      just        = c("center", "top"),
      default.units= "inches"
    )
    
    # Title and legend
    plotText(
      label         = paste("QTL:", trait, "(", norm, "/", drug, ") -",
                            chr_str, "Peak:", round(peak_pos, 2), "Mb"),
      x             = PLOT_PARAMS$x,
      y             = 0.05,
      just          = c("center", "top"),
      fontface      = "bold",
      fontsize      = 12,
      default.units = "inches"
    )
    plotLegend(
      legend        = c("Original LOD", "pyResults"),
      fill          = c(COLORS[[norm]]["Sig"], "black"),
      border        = FALSE,
      x             = PLOT_PARAMS$x,
      y             = 0,
      width         = 1,
      height        = 0.75,
      just          = c("center", "top"),
      fontsize      = 10,
      default.units = "inches",
      orientation   = "h"
    )
    
    dev.off()
  })

print("Finished plotting significant loci with pyResults sub-plots.")
