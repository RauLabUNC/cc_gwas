# Load necessary libraries
library(tidyverse)
library(plotgardener)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(miqtl)

merge_blocks <- function(blocks) {
  merged <- blocks %>%
    arrange(chr, lower.pos, upper.pos) %>%
    mutate(group = cumsum(!(chr == lag(chr, default = data.table::first(chr)) & 
                              (lower.pos == lag(lower.pos, default = data.table::first(lower.pos)) | 
                                 upper.pos == lag(upper.pos, default = data.table::first(upper.pos)))))) %>%
    group_by(group) %>%
    summarise(
      chr = data.table::first(chr),
      lower.pos = min(lower.pos),
      upper.pos = max(upper.pos),
      block = min(block)
    ) %>%
    ungroup() %>%
    mutate(block = row_number()) %>%
    dplyr::select(block, chr, lower.pos, upper.pos)
  return(merged)
}

# Set arguments
args <- c("CSA", "iso")
#args <- commandArgs(trailingOnly = TRUE)
#phenotype_of_interest <- args[1]

# Load genome assembly
mm39 <- assembly(Genome = "mm39_GRCm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene", OrgDb = "org.Mm.eg.db")

# Read and process threshold data
threshold_file <- file.path("data/processed/scan_thresholds", args[2], paste0(args[1], "_scan.rds"))
threshold <- readRDS(threshold_file)
threshold <- get.gev.thresholds(threshold, use.lod = T, percentile = 0.85, type = "min")

# Read and process range data
range_file <- file.path("data/processed/sig_loci/oneFiveLODRanges", args[2], paste0(args[1], ".csv"))
ranges <- read.csv(range_file)

# Define the groups to plot
blocks <- ranges %>%
  filter(!is.na(block)) %>%
  group_by(block) %>%
  slice_head() %>%
  dplyr::select(block, chr, lower.pos, upper.pos) %>%
  ungroup() |> as.data.frame()

if(nrow(blocks) != 0){
# Merge blocks
merged_blocks <- merge_blocks(blocks) |> as.data.frame()

# Select a block to plot

for(i in unique(merged_blocks$block)) {
  
block <- merged_blocks %>% filter(block == i)
block.chr <- block$chr
range <- ranges %>% filter(chr == block.chr)

# Convert pos from Mb to bp
block.chr <- paste0("chr", block.chr)
range <- range %>%
  mutate(pos = pos * 1e6,
         chrom = paste0("chr", chr),
         p = 10^(-LOD)) %>%
  dplyr::select(chrom, pos, p)

block <- block %>%
  mutate(lower.pos = lower.pos * 1e6,
         upper.pos = upper.pos * 1e6)

# Set x axis range based on lower/upper.pos
delta.pos <- block$lower.pos - block$upper.pos
upper.lim <- as.integer(block$upper.pos - 0.8 * delta.pos)
lower.lim <- as.integer(block$lower.pos + 0.8 * delta.pos)
# Create the LOD score plot

# Plot the results
output_dir <- file.path("results/locusZoomPlots", args[[2]], args[[1]])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

png_name <- paste0(output_dir, "/block_", i, ".png")
png(file = png_name,
    width = 7, 
    height = 6,
    units = "in",
    res = 600)


# Create a plotgardener page
pageCreate(
  width = 6, height = 5, default.units = "inches",
  showGuides = FALSE, xgrid = 0, ygrid = 0
)

params <- pgParams(assembly = mm39, quiet = T, 
                   just = c("center", "center"), 
                   default.units = "inches",
                   chromstart = upper.lim, chromend = lower.lim)
# Plot LOD score data using plotSignal
manhattanPlot <- plotManhattan(
  data = range,
  chrom = block.chr, 
  trans = "-log10", sigVal = 10^(-threshold),
  x = 3, y = 1, width = 6, height = 2,
  just = c("center", "top"), 
  xfield = "pos", yfield = "p", fill = "black", sigCol = "darkorange",
  sigLine = T, baseline = T,
  yscale_reverse = F, params = params
)

## Annotate y-axis
annoYaxis(
  plot = manhattanPlot, at = seq(0, max(ceiling(-log10(range$p)))),
  axisLine = TRUE, fontsize = 8, params = params)

annoHighlight(
  plot = manhattanPlot,
  chrom = block.chr,
  chromstart = block$upper.pos, chromend = block$lower.pos,
  x = 3, y = 1, width = 6, height = 2,
  just = c("center", "top"), 
  params = params
)

## Plot text label
## Plot text, adjusting fontsize and fontface
plotText(
  label = "1.5 LOD interval range", fontsize = 12, 
  x = 3, y = 0.9, just = c("center", "center"), default.units = "inches", color = "grey",
  params = params
)

plotText(
  label = paste0("Block #", i, " of ", nrow(merged_blocks)), fontsize = 12, fontface = "bold",
  x = -0.3, y =0.8, just = c("left", "center"), default.units = "inches", color = "black",
  params = params
)

plotText(
  label = paste0("QTL for ", args[[1]], " in ", args[[2]], " CC mice"), fontsize = 14, fontface = "bold",
  x = 3, y = 0.5, just = c("center", "center"), default.units = "inches", color = "black",
  params = params
)

plotGenes(
  chrom = block.chr, chromstart = upper.lim, chromend = lower.lim,
  assembly = mm39,
  x = 3, y = 3, width = 6, height = 2,
  just = c("center", "top"), default.units = "inches",
  params = params
)

plotGenomeLabel(
  chrom = block.chr, chromstart = upper.lim, chromend = lower.lim,
  assembly = mm39,
  x = 3, y = 5, length = 6, scale = "Mb",
  just = c("center", "top"), default.units = "inches",
  params = params, quiet = T
)
dev.off()

}
}

# Temp solution to issue of uncertain script outputs in snakemake
output_dir <- file.path("results/locusZoomPlots", args[[2]], args[[1]])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.table("dummy", paste0(output_dir, "/", "dummy.txt"))
