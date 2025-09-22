suppressPackageStartupMessages({
  library(tidyverse)
  library(plotgardener)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(GenomeInfoDb)
  library(AnnotationHub)
  library(RColorBrewer)
})


mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

# Load in annotated loci to groups
loci <- read.csv("results/sig_regions/locus_clusters.csv") |> as_tibble()
loci_grouped <- loci %>%
  group_by(locus_cluster) %>%
  summarise(
    chr = dplyr::first(chr),
    chr_formatted = paste0("chr", dplyr::first(chr)),
    start_bp = min(as.numeric(start_bp)),
    end_bp = max(as.numeric(end_bp)),
    trait_drug = paste(paste(trait, drug, sep = ":"), collapse = "/"),
    type = ifelse(n() == 1, "Monotropic", "Pleiotropic"),
    .groups = "drop"
  )
# Desired chromosome order: 1–19, X, Y
target_order <- c(as.character(1:19), "X", "Y")

# Normalize loci$chr to uppercase character (no "chr" prefix)
loc_chr <- toupper(as.character(loci$chr))

# Determine which targets are present in loci and in the genome
tx_db <- TxDb.Mmusculus.UCSC.mm39.knownGene
chromSizes <- GenomeInfoDb::seqlengths(tx_db)
maxChromSize <- max(chromSizes)

ordered_chr <- paste0("chr", target_order)
is_present <- target_order %in% loc_chr

# Count contiguous missing-chromosome groups
missing_rle <- rle(!is_present)
n_missing_groups <- if (length(missing_rle$values)) sum(missing_rle$values) else 0

# Report which targets are missing (will be plotted as small gaps)
missing_chrs <- target_order[!is_present]
if (length(missing_chrs) > 0) {
  message("Missing chromosomes (plotting as gaps): ", paste(missing_chrs, collapse = ", "))
}

# set up spacing and dimensions
width <- 3.4
height <- 2

# Adjustable gaps (in inches)
inter_chr_gap <- 0.05    # space between plotted chromosomes
missing_chr_gap <- 0.05  # space for each missing chromosome

# Solve for chrom_width so layout fits exactly into 'width'
n_present <- sum(is_present)
if (n_present == 0) stop("No chromosomes from target_order present in loci and genome.")

last_is_present <- tail(is_present, 1)
B <- (n_present * inter_chr_gap) + (n_missing_groups * missing_chr_gap) - ifelse(last_is_present, inter_chr_gap, 0)
chrom_width <- (width - B) / n_present
if (chrom_width <= 0) stop("Layout over-constrained: decrease gaps or increase width.")

# Vertical scale respects relative chromosome lengths
label_y_offset <- 0.1
usable_height <- height - label_y_offset

pdf("results/paper_figures/fig2b.pdf", width = height + 0.2, height = width + 0.1, useDingbats = FALSE)
pageCreate(
  width = height + 0.2, height = width + 0.1, default.units = "inches",
  showGuides = FALSE, xgrid = 0, ygrid = 0
)
# March across the target_order, plotting present chromosomes and
# advancing x by a small gap for missing ones.
x_current <- 0.05
for (i in seq_along(target_order)) {
  #i <- 1
  chr_label <- ordered_chr[i]
  if (is_present[i]) {
    ideo_height <- (usable_height * chromSizes[[chr_label]]) / maxChromSize
    x_center <- x_current + (chrom_width / 2)

    ideo_plot <- plotIdeogram(
      chrom = chr_label, assembly = "mm10",
      showBands = FALSE,
      orientation = "h",
      x = 0.2, y = x_center,
      width = ideo_height, height = chrom_width,
      just = "left"
    )
    # Loop to add all monotropic and pleiotropic loci on this chromosome
    loci_this_chr <- loci_grouped %>% filter(chr_formatted == chr_label)
    if (nrow(loci_this_chr) > 0) {
      for (j in seq_len(nrow(loci_this_chr))) {
        locus <- loci_this_chr[j, ]
        region <- pgParams(chrom = locus$chr_formatted, 
                           chromstart = locus$start_bp, 
                           chromend = locus$end_bp)
        if(locus$type == "Pleiotropic") {
          fill_color <- "#1B9E77FF"
        } else {
          fill_color <- "#D95F02FF"
        }
        annoHighlight(
          plot = ideo_plot, params = region,
          fill = fill_color,
          y = x_center, height = chrom_width, just = c("center"), default.units = "inches"
        )
      }
    }
    plotText(
      label = gsub("chr", "", chr_label),
      x = 0.05, y = x_center, fontsize = 8, rot = 270,
      just = c("bottom", "center")
    )
    # advance by chromosome width + inter-chromosome gap
    x_current <- x_current + chrom_width + inter_chr_gap
  } else {
    # advance once per contiguous block of missing chromosomes
    if (i == 1 || is_present[i - 1]) {
      x_current <- x_current + missing_chr_gap
    }
  }
}
dev.off()





