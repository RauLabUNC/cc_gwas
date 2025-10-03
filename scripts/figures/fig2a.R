suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(plotgardener)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(miqtl)
  library(igraph)
  library(GenomicRanges)
  library(RColorBrewer)
  library(ragg)
  library(png)
})

# --- Data ---
scan <- readRDS("data/processed/ropscan/TH_Ctrl.rds")
threshold <- readRDS("data/processed/scan_thresholds/TH_Ctrl_threshold.rds")

miqtl_df_for_plot <- tibble(
  snp = names(scan$LOD),
  chr = as.character(scan$chr),
  pos = scan$pos$Mb * 1e6,
  lod = scan$LOD
) |>
  filter(!is.na(pos), !is.na(lod)) |>
  # plotgardener expects a 'p' to transform; we feed 10^-LOD so -log10(p) == LOD
  transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod), snp)

# --- Assembly ---
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

# --- Figure sizing & layout (no outer margins) ---
PLOT_DIMS <- list(page_width = 3.5, page_height = 2)  # inches
# fill the page completely: x=0, y=0, width=page_width, height=page_height
PLOT_PARAMS <- list(x = 0.32, y = 0.15, plot_width = PLOT_DIMS$page_width-0.32, plot_height = PLOT_DIMS$page_height-0.3)

# Alternating greys by chromosome
chr_fill <- c("grey85", "grey60")

# Y range in LOD units (since -log10(p) == LOD for p = 10^-LOD)
miqtl_ylim <- c(
  0,
  max(c(-log10(miqtl_df_for_plot$p), threshold, 5), na.rm = TRUE)
)
# 1) RASTERIZE THE DOTS (and baseline/sig line) ONLY  -------------------------
dots_png <- "results/paper_figures/fig2a_points.png"
# --- Output ---
if (!dir.exists("results/paper_figures")) dir.create("results/paper_figures", recursive = TRUE)
ragg::agg_png(
  dots_png,
  width = PLOT_DIMS$page_width,
  height = PLOT_DIMS$page_height,
  units = "in", res = 1200, background = "transparent"
)

# Create a full-bleed page (no extra space)
pageCreate(
  width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height,
  default.units = "inches", showGuides = FALSE
)
# --- Manhattan plot ---
miqtl_plot <- plotManhattan(
  data = miqtl_df_for_plot,
  assembly = mm39,
  range = miqtl_ylim,
  trans = "-log10",
  sigVal = 10^(-threshold),     # threshold converted to "p"
  fill   = chr_fill,            # alternate greys per chromosome
  sigCol = "#1f78b4",           # significant points drawn in blue
  sigLine = TRUE,
  baseline = TRUE,
  cex = 0.15,                   # tuned for 3.4" x 2.5"
  space = 0.02,                 # minimal spacing between chromosomes
  x = PLOT_PARAMS$x, y = PLOT_PARAMS$y,
  width = PLOT_PARAMS$plot_width, height = PLOT_PARAMS$plot_height,
  just = c("left", "top"),
  default.units = "inches",
  draw = TRUE
)
dev.off()
raster_plot <- readPNG("results/paper_figures/fig2a_points.png")

# 2) BUILD THE FINAL PDF WITH VECTOR AXES/LABELS OVER THE RASTER --------------
pdf("results/paper_figures/fig2a.pdf", width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height, useDingbats = FALSE)
pageCreate(width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height, default.units = "inches", showGuides = FALSE)

# Place the rasterized dots to fill the page
plotRaster(
  image = raster_plot,
  x = 0, y = 0,
  width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height,
  just = c("left","top"),
  default.units = "inches"
)
man_dummy <- plotManhattan(
  data = miqtl_df_for_plot, assembly = mm39,
  range = miqtl_ylim, trans = "-log10",
  x = PLOT_PARAMS$x, y = PLOT_PARAMS$y,
  width = PLOT_PARAMS$plot_width, height = PLOT_PARAMS$plot_height,
  just = c("left","top"),
  draw = FALSE
)
# --- Y axis (8 pt) drawn inside the plot area ---
annoYaxis(
  plot = man_dummy,
  at = pretty(miqtl_ylim),
  axisLine = TRUE,
  fontsize = 8,
  main = TRUE
)

# --- Y label (8 pt), placed inside left side to avoid outer margins ---
plotText(
  label = "LOD",
  x = 0,                                   # slight inset so it doesn't clip
  y = PLOT_PARAMS$y + PLOT_PARAMS$plot_height / 2,
  rot = 90,
  fontsize = 8,
  just = c("center", "top"),
  default.units = "inches"
)

# --- Chromosome labels at the top edge (8 pt), inside the plotting area ---
annoGenomeLabel(
  plot = man_dummy,
  x = PLOT_PARAMS$x,
  y = PLOT_PARAMS$y + PLOT_PARAMS$plot_height,
  fontsize = 7,
  just = c("left", "top"),
  default.units = "inches"
)

dev.off()







