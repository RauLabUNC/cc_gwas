suppressPackageStartupMessages({
  library(tidyverse)
  library(plotgardener)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(GenomeInfoDb)
  library(AnnotationHub)
  library(RColorBrewer)
  library(data.table)
})


#### Load inputs for what we're plotting
# Scan and threshold for Chr 10 LV.Vold.21 ISO
scan <- readRDS("data/processed/ropscan/LV.Vold.21_Iso.rds")
alt_scan <- readRDS("data/processed/ropscan/LV.Vold.21_Ctrl.rds")
threshold <- readRDS("data/processed/scan_thresholds/LV.Vold.21_Iso_threshold.rds")
alt_threshold <- readRDS("data/processed/scan_thresholds/LV.Vold.21_Ctrl_threshold.rds")

# List of genes
genes <- read.csv("data/processed/joinLoci/relational_tables/genes_mouse.csv")

# Loci information
loci <- read.csv("results/sig_regions/locus_clusters.csv") |> as_tibble()

#### Data wrangling 
miqtl_df_for_plot <- tibble(
  snp = names(scan$LOD),
  chr = as.character(scan$chr),
  pos = scan$pos$Mb * 1e6,
  lod = scan$LOD
) |>
  filter(!is.na(pos), !is.na(lod)) |>
  # plotgardener expects a 'p' to transform; we feed 10^-LOD so -log10(p) == LOD
  transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod), snp)

alt_miqtl_df_for_plot <- tibble(
  snp = names(alt_scan$LOD),
  chr = as.character(alt_scan$chr),
  pos = alt_scan$pos$Mb * 1e6,
  lod = alt_scan$LOD
) |>
  filter(!is.na(pos), !is.na(lod)) |>
  transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod), snp)

# Find loci that overlap with our region of interest
target_cluster <- loci %>%
  filter(trait == "LV.Vold.21" & drug == "Iso") %>%
  dplyr::pull(locus_cluster) %>%
  unique()
overlaps <- loci |> 
    filter(locus_cluster == target_cluster) |>
    mutate(chrom = paste0("chr", chr),
           start = start_bp,
           end   = end_bp, 
           strand= "-",
           traitXdrug = paste0(trait, ": ", drug)) |>  
    dplyr::select(chrom, start, end, strand, trait, drug, traitXdrug)

#### PlotGardener params and objects
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

miqtl_ylim <- c(0, max(c(-log10(miqtl_df_for_plot$p), threshold, 5), na.rm = TRUE))
plot_start_bp <- 27000000 #max(overlaps$start)
plot_end_bp <- 36000000 #min(overlaps$end) 
bounds_bp <- c(plot_start_bp, plot_end_bp)
plot_chr <- unique(overlaps$chrom)

PLOT_DIMS <- list(page_width = 3.5, page_height = 3.5)  # inches
PLOT_PARAMS <- list(x = 0.35, plot_width = PLOT_DIMS$page_width-0.4,
                    y1 = 0.2, height1 = 1.2, # Manhattan
                    y2 = 1.6, height2 = 1.2) # Founder effects, y2 = y1 + height1 + 0.1


# --- plotGardener setup ---
pdf("results/paper_figures/fig2c.pdf", width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height)
pageCreate(width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height, default.units = "inches", showGuides = FALSE)

params_genome <- pgParams(
assembly   = mm39, 
chrom      = plot_chr,
chromstart = plot_start_bp,
chromend   = plot_end_bp
)

## Manhattan plots ##
alt_miqtl_plot <- plotManhattan(
    data = alt_miqtl_df_for_plot, params = params_genome,
    range = miqtl_ylim, trans = "-log10", sigVal = 10^(-alt_threshold),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$y1 ,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$height1,
    just = c("left", "top"), xfield = "", yfield = "",
    fill = "#e8b969ff", sigCol = "#b4761fff", sigLine = FALSE, baseline = FALSE,
    default.units = "inches"
    )
miqtl_plot <- plotManhattan(
    data = miqtl_df_for_plot, params = params_genome,
    range = miqtl_ylim, trans = "-log10", sigVal = 10^(-threshold),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$y1 ,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$height1,
    just = c("left", "top"), xfield = "pos", yfield = "p",
    fill = "#a6cee3", sigCol = "#1f78b4", sigLine = TRUE, baseline = TRUE,
    default.units = "inches"
    )
annoYaxis(plot = miqtl_plot, at = pretty(miqtl_ylim),
        axisLine     = TRUE,
        fontsize     = 8,
        main         = TRUE)
plotText(label = "LOD",  
        x       = 0.05,
        y       = PLOT_PARAMS$y1 + PLOT_PARAMS$height1 / 2,
        rot     = 90,
        fontsize= 8,
        just    = c("center", "top"),
        default.units = "in")

## Founder effects plot ##
pos_in_range_logical <- scan$pos$Mb * 10^6 >= (plot_start_bp - 10^5) & scan$pos$Mb * 10^6 <= (plot_end_bp + 10^5)   
chr_logical <- scan$chr == gsub("chr", "", plot_chr)
marker_positions_bp <- scan$pos$Mb[pos_in_range_logical & chr_logical] * 10^6
  
  markers <- data.frame(marker = scan$loci[pos_in_range_logical & chr_logical],
                        start  = marker_positions_bp) |> 
    arrange(start)
  
  allele_effects_matrix <- scan$allele.effects[,markers$marker, drop = FALSE] # drop=FALSE to handle single marker case
  allele_effects_transposed <- t(scan$allele.effects) |> as.data.frame()
  allele_effects_transposed$marker <- rownames(allele_effects_transposed)
  
  
  num_strains <- ncol(allele_effects_transposed) -1 
  founder_strains <- rownames(scan$allele.effects)
  
plot_data_list <- lapply(1:num_strains, function(strain){
    curr_strain <- colnames(allele_effects_transposed)[[strain]]
    temp_df <- allele_effects_transposed |> 
        dplyr::select(marker, founder_strains[strain]) %>% 
        filter(marker %in% markers$marker)

    temp_df <- left_join(temp_df, markers, by = "marker")

    temp_df <- temp_df |> 
        mutate(chrom = plot_chr) |> 
        arrange(start)
    colnames(temp_df)[2] <- "score"

    # Need to change data points into ranges, so make a p value continue downstream until the next 
    temp_df$end <- c(temp_df$start[2:nrow(temp_df)] - 1L, 
                        temp_df$start[  nrow(temp_df)] + 1)   
    temp_df <- temp_df |> 
        dplyr::select(chrom, start, end, score)
    return(temp_df)
    }
)
names(plot_data_list) <- founder_strains


# Get the max abs haplotype effects to set y axis 
pos_allele_max <-  plot_data_list |> rbindlist() |> pull(score) |> max() |> ceiling()
neg_allele_min <-  plot_data_list |> rbindlist() |> pull(score) |> min() |> floor()

strain_colors <- rep(c("#ff0", "#888", "#f88", "#11f", "#0cf", "#0a0", "#f00", "#90e"), length.out = num_strains)
founders <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/H1LtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")

signalPlots <- c()

#------- Plot the strain effects ------#
signalPlots[[1]] <- plotgardener::plotSignal(
    data = plot_data_list[[1]],
    params = params_genome,
    range = c(neg_allele_min, pos_allele_max), 
    linecolor = strain_colors[[1]],
    fill = NA, 
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$y2,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$height2,
    just = c("left", "top"), default.units = "inches",
    baseline = TRUE, 
    baseline.color = "grey")

annoYaxis(plot = signalPlots[[1]], 
        at = c(neg_allele_min, 0, pos_allele_max),
        axisLine     = TRUE,
        fontsize     = 8,
        main         = TRUE)

lapply(2:8, function(strain_data){
signalPlots[[strain_data]] <- plotgardener::plotSignal( 
    data = plot_data_list[[strain_data]],
    params = params_genome,
    range = c(neg_allele_min, pos_allele_max), 
    linecolor = strain_colors[strain_data],
    fill = NA, 
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$y2,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$height2,
    just = c("left", "top"), default.units = "inches",
    baseline = FALSE)
}
)

plotText(label = "Founder Effects",  
        x       = 0.05,
        y       = PLOT_PARAMS$y2 + PLOT_PARAMS$height2 / 2,
        rot     = 90,
        fontsize= 8,
        just    = c("center","top"),
        default.units = "in")

# Get order of strains by effect at the end of the locus
strain_end_order <- lapply(names(plot_data_list), function(strain){
    df <- plot_data_list[[strain]] |> filter(start < plot_end_bp & end > plot_end_bp) |> slice_head(n=1) 
    return(df)
})


strain_end_order_df <- rbindlist(strain_end_order)
strain_end_order_df$strain <- founders

strain_end_order <- strain_end_order_df |> arrange(-score) |> dplyr::pull(strain)

# Reorder strain_colors to match the new strain order
strain_colors <- strain_colors[match(strain_end_order, founders)]

  legendPlot <- plotLegend(
    legend = strain_end_order[1:4],
    fill = strain_colors[1:4],
    border = FALSE, 
    x = 0.6, 
    y = PLOT_DIMS$page_height ,
    width = 0.5, 
    height = 0.5,
    fontsize = 8,
    just = c("left", "bottom"),
    orientation = "v",
    default.units = "inches"
  )
  legendPlot <- plotLegend(
    legend = strain_end_order[5:8],
    fill = strain_colors[5:8],
    border = FALSE, 
    x = 2, 
    y = PLOT_DIMS$page_height ,
    width = 0.5, 
    height = 0.5,
    fontsize = 8,
    just = c("left", "bottom"),
    orientation = "v",
    default.units = "inches"
  )

  legendPlot <- plotLegend(
    legend = c("Iso", "Ctrl"),
    fill = c("#1f78b4", "#b4761fff"),
    border = FALSE, 
    x = PLOT_PARAMS$x + PLOT_PARAMS$plot_width/2, 
    y = 0.05, 
    width = 0.5, 
    height = 0.1,
    fontsize = 8,
    just = c("center", "top"),
    orientation = "h",
    default.units = "inches"
  )

  
  # --- Add Genome Label ---
  annoGenomeLabel(
    plot = signalPlots[[1]], 
    params = params_genome,
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$y2 + PLOT_PARAMS$height2,
    scale = "Mb", fontsize = 8,
    just = c("left", "top"), default.units = "inches")

dev.off()




# --- plotGardener setup ---
pdf("results/paper_figures/fig2b_scalebar.pdf", width = 3.5, height = 2.2)
pageCreate(width = 3.5, height = 2.2, default.units = "inches", showGuides = FALSE)

maxChromSize <- 195154279 # from fig2b.R
ideo_height <- 1.9

# Make a signal plot with a max value of 19.5 Mb, 1.9 tall
dummy_plot <- plotgardener::plotSignal(
    data = plot_data_list[[1]],
    params = params_genome,
    range = c(0, 19.5), 
    linecolor = strain_colors[[1]],
    fill = NA, 
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$y1,
    width          = PLOT_PARAMS$plot_width,
    height         = 1.9,
    just = c("left", "top"), default.units = "inches",
    baseline = TRUE, 
    baseline.color = "grey", draw = FALSE)
annoYaxis(plot = dummy_plot, 
        at = c(0, 10, 20),
        axisLine     = TRUE,
        fontsize     = 8,
        main         = TRUE)
plotText(label = "Megabases",  
        x       = 0,
        y       = PLOT_PARAMS$y1 + 1.9 / 2,
        rot     = 90,
        fontsize= 8,
        just    = c("center","top"),
        default.units = "in")
dev.off()

