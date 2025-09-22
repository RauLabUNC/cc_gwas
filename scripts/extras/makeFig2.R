suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(plotgardener) # For locus plots
  library(TxDb.Mmusculus.UCSC.mm39.knownGene) 
  library(org.Mm.eg.db) # For gene ID mapping
  library(openxlsx) # For writing Excel files
  library(miqtl)
  library(igraph)         
  library(GenomicRanges)
  library(RColorBrewer)
})

# Mimic optparse behavior
opt <- list()

opt$input_sig_regions <- "data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv"
opt$input_all_scans <- "data/processed/trait_qtl/all_scans.rds"
opt$input_all_thresholds <- "data/processed/trait_qtl/all_thresholds.rds"
opt$input_merged_gene_info <- "data/processed/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv"

base_relational <- "data/processed/joinLoci/relational_tables"
paths_multitrait <- list(
  genes_mouse   = fs::path(base_relational, "genes_mouse.csv"),
  orthology     = fs::path(base_relational, "orthology.csv"),
  associations  = fs::path(base_relational, "associations.csv"),
  trait_loci    = fs::path(base_relational, "traitLoci.csv"),
  mouse_pheno   = fs::path(base_relational, "mouseGenePhenotypes.csv")
)

genes_mouse   <- fread(paths_multitrait$genes_mouse)
#ortho_mouse2h <- fread(paths_multitrait$orthology)
trait_loci    <- fread(paths_multitrait$trait_loci) # Contains miQTL trait loci
#associations  <- fread(paths_multitrait$associations) # Human disease associations
#mouse_pheno   <- fread(paths_multitrait$mouse_pheno)

# === Data used in plotLoci ===
sig_regions    <- read_csv(opt$input_sig_regions) %>%
  mutate(
    across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
    chr = as.character(chr)
  ) %>%
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug)

scan_data      <- readRDS(opt$input_all_scans) # miQTL scan results
threshold_data <- readRDS(opt$input_all_thresholds)

# === Merged table output from multTrait_cis-eQTL_nrvmExp.R ===
merged_gene_info <- fread(opt$input_merged_gene_info)

# === CC Founder Variant Information ===
muNoSplice <- fread("data/processed/joinLoci/relational_tables/ccVariants/Gene_Mutations_CC_WT_NoSplice_OnlyDeleterious_3_31_25.csv")
muDriving  <- fread("data/processed/joinLoci/relational_tables/ccVariants/Driving_Mutations_by_SNP_Filtered.csv")

# === Bulk RNAseq, CPM and VST, and DESeq2 ===
#rna_info <- fread("data/processed/joinLoci/bulk_exp/5d_VST_Info_250429.csv", drop = 1)
#cpm_ctrl <- read.csv("data/processed/joinLoci/bulk_exp/Ctrl_CPM.csv")
#cpm_iso  <- read.csv("data/processed/joinLoci/bulk_exp/Iso_CPM.csv")

#deseq <- read.csv("data/processed/joinLoci/bulk_exp/deseq2_13kGenes_SexCovar_250521.csv", row.names = 1)

# === eQTL mapping results ===
cis_eQTL <- read.csv("data/processed/joinLoci/bulk_exp/VST_1or5Mb_Cis_eQTL_simplified.csv", header = T, row.names = 1)

# === get the genome assembly

# Genome assembly
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

message("All data loaded.")

### Define plotting parameters ###
params <- c()
# Constants
COLORS <- c(Sig = "#1f78b4", nonSig = "#a6cee3")
PLOT_DIMS <- list(page_width = 6, page_height = 3, res = 300)
PLOT_PARAMS <- list(x = 3, plot_width = 5, plot_height = 1.5, plot_y = 0.5)
GENE_DIMS <- list(height = 2, y_offset = 0.5, label_offset = 0.1)
MIN_YLIM <- 5
#strain_colors <- rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), length.out = num_strains)
founders <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/H1LtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")

# Identify locus groups #####
sig_regions    <- read_csv(opt$input_sig_regions) |>
  mutate( across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
        chr = as.character(chr),
        start_bp = as.integer(floor(upper_pos_lod_drop * 1e6)),
        end_bp   = as.integer(floor(lower_pos_lod_drop * 1e6))) |>
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug) |>
  dplyr::select(-upper_pos_lod_drop, -lower_pos_lod_drop) |>
  arrange(chr, start_bp, end_bp)

sig_regions[1,]
seqs <- ifelse(grepl("^chr", sig_regions$chr), sig_regions$chr, paste0("chr", sig_regions$chr))
gr <- GRanges(seqnames = seqs,
              ranges   = IRanges(sig_regions$start_bp, sig_regions$end_bp),
              strand   = "*",
              mcols    = sig_regions)
str(gr)
hits <- GenomicRanges::findOverlaps(gr) 
hits
edges <- as.data.frame(hits)[ , 1:2]               # queryHits | subjectHits
edges
edges <- subset(edges, queryHits != subjectHits)   # drop the trivial self‑edges
edges <- unique(edges)                             # remove duplicated pairs

g <- graph_from_data_frame(edges,
                           directed = FALSE,
                           vertices = data.frame(name = seq_along(gr)))

## Connected components = overlapping loci 
comp <- components(g)$membership
sig_regions$locus_cluster <- comp

viewable <- sig_regions |> arrange(locus_cluster, chr, peak_pos) |> dplyr::select(locus_cluster, chr, start_bp, end_bp, trait, drug, peak_pos, max_lod)
print(viewable, n=50)
write.csv(viewable, "results/sig_regions/locus_clusters.csv", row.names = FALSE)

# unique trait-drug combos, then keep traits that appear with exactly one drug
trait_specific <- viewable %>%
    dplyr::distinct(trait, drug) %>%
    dplyr::add_count(trait, name = "n_drugs") %>%
    dplyr::filter(n_drugs == 1) %>%
    dplyr::select(trait, drug) %>%
    dplyr::filter(!grepl("\\.0", trait))

multi_locus <- viewable %>%
  dplyr::distinct(trait, locus_cluster) %>%
  dplyr::add_count(trait, name = "n_clust") %>%
  dplyr::filter(n_clust > 1) %>%
  dplyr::select(trait, locus_cluster) %>%
  dplyr::filter(!grepl("\\.0", trait))

multi_chr <- viewable %>%
  dplyr::distinct(trait, chr) %>%
  dplyr::add_count(trait, name = "n_chr") %>%
  dplyr::filter(n_chr > 1) %>%
  dplyr::select(trait, chr) %>%
  dplyr::filter(!grepl("\\.0", trait))

# Get loci that appear in mulit-trait loci
# Count loci per cluster
cluster_sizes <- viewable %>%
  dplyr::filter(!is.na(locus_cluster)) %>%
  dplyr::count(locus_cluster, name = "n_loci")

# How many clusters have >= 2 loci
n_multi_clusters <- cluster_sizes %>% dplyr::filter(n_loci >= 2) %>% nrow()
mono_clusters <- cluster_sizes %>% dplyr::filter(n_loci == 1) %>% nrow()
# How many individual loci are in those clusters
total_loci_in_multi <- cluster_sizes %>% dplyr::filter(n_loci >= 2) %>% dplyr::summarise(sum(n_loci)) %>% dplyr::pull()
###########################################33
# --- 5. Main Loop to Generate Packets --- ####
##########################################
# Define which loci to process.
group_to_process <- sig_regions |>
    dplyr::filter(locus_cluster == 4)

# Make the overlaps dataframe to plot ranges of loci within loci cluster
overlaps <- group_to_process |> 
    mutate(chrom = paste0("chr", chr),
           start = start_bp,
           end   = end_bp, 
           strand= "-",
           traitXdrug = paste0(trait, ": ", drug)) |>  
    dplyr::select(chrom, start, end, strand, trait, drug, traitXdrug)

# locus_info should contain: chr, peak_pos, lower_pos_lod_drop, upper_pos_lod_drop, trait, drug
output_dir <- "results/sig_regions/fig2"

current_chr <- group_to_process$chr[1]
current_start <- min(group_to_process$start_bp)
current_end   <- max(group_to_process$end_bp)

genes_in_locus <- genes_mouse |>
  dplyr::filter(chr == current_chr &
                start_bp >= current_start &
                end_bp <= current_end) %>%
  dplyr::arrange(start_bp)
  
  # Ensure positions are in base pairs for plotgardener
  plot_start_bp <- current_start - 5*1e5 # Add padding
  plot_end_bp   <- current_end + 5*1e5 # Add padding
  bounds_bp <- c(plot_start_bp, plot_end_bp)
  
  # Ensure start is not negative
  plot_start_bp <- max(0, plot_start_bp)
  
  dir_path <- file.path(output_dir, "zoomPlots")
  plot_file_name <- file.path(dir_path, "chr10_2.pdf")
  if(!dir.exists(dir_path)){dir.create(dir_path)}
  message("Generating locus zoom plot: ", plot_file_name)
  
  ######################################
  #### --- plotGardener setup --- #####
  ######################################
  pdf(plot_file_name, width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height)
  pageCreate(width = PLOT_DIMS$page_width, height = PLOT_DIMS$page_height, default.units = "inches", showGuides = FALSE)

  params_genome <- pgParams(
    assembly   = mm39#, 
    chrom      = paste0("chr", current_chr),
    chromstart = plot_start_bp,
    chromend   = plot_end_bp
  )

  # stand in for a loop through all loci in group
  locus_info <- group_to_process[2,]

  # --- Plot miQTL LOD scores ---
  current_scan_data <- scan_data[[paste0(locus_info$trait, "_", locus_info$drug)]]
  
  if(locus_info$drug == "Ctrl"){
    alt_drug <- "Iso"
  } else if(locus_info$drug == "Iso"){
    alt_drug <- "Ctrl"
  }
  alt_scan_data <- scan_data[[paste0(locus_info$trait, "_", alt_drug)]]


  miqtl_df_for_plot <- tibble(
    marker = names(current_scan_data$LOD),
    chr    = as.character(current_scan_data$chr),
    pos    = current_scan_data$pos$Mb * 1e6, 
    lod    = current_scan_data$LOD
  ) %>%
    filter(chr == current_chr, !is.na(pos), !is.na(lod)) %>% # Filter for current chromosome
    transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod)) #p doesn't really mean P here, plotGardener just wants to take the log10 of this number
  alt_miqtl_df_for_plot <- tibble(
    marker = names(alt_scan_data$LOD),
    chr    = as.character(alt_scan_data$chr),
    pos    = alt_scan_data$pos$Mb * 1e6, 
    lod    = alt_scan_data$LOD
  ) %>%
    filter(chr == current_chr, !is.na(pos), !is.na(lod)) %>% # Filter for current chromosome
    transmute(chrom = paste0("chr", chr), pos, p = 10^(-lod)) #p doesn't really mean P here, plotGardener just wants to take the log10 of this number
  
  # Determine y-axis limits for miQTL plot
  alt_miqtl_threshold_val <- threshold_data[[paste0(locus_info$trait, "_", alt_drug, "_threshold")]]
  miqtl_threshold_val <- threshold_data[[paste0(locus_info$trait, "_", locus_info$drug, "_threshold")]]
  miqtl_ylim <- c(0, max(c(-log10(miqtl_df_for_plot$p), miqtl_threshold_val, 5), na.rm = TRUE) + 1)
  
  ## Plot miQTL ##
  alt_miqtl_plot <- plotManhattan(
    data = alt_miqtl_df_for_plot, params = params_genome,
    range = miqtl_ylim, trans = "-log10", sigVal = 10^(-alt_miqtl_threshold_val),
    x              = PLOT_PARAMS$x,
    y              = PLOT_PARAMS$plot_y ,
    width          = PLOT_PARAMS$plot_width,
    height         = PLOT_PARAMS$plot_height,
    just = c("center", "top"), xfield = "", yfield = "",
    fill = "#e8b969ff", sigCol = "#b4761fff", sigLine = FALSE, baseline = FALSE,
    default.units = "inches"
  )
  
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
  plotText(label = "LOD",  
           x       = PLOT_PARAMS$plot_width/2 + PLOT_PARAMS$x + 0.4,
           y       = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height / 2,
           rot     = 270,
           fontsize= 8,
           just    = "center",
           default.units = "in")
  legendPlot <- plotLegend(
    legend = c("Iso", "Control"),
    fill = c("#a6cee3", "#e8b969ff"),
    border = FALSE, 
    x = PLOT_PARAMS$x, y = 0.5, width = 1.2, height = 0.5,
    fontsize = 8,
    just = c("center", "center"),
    orientation = "h",
    default.units = "inches"
  )

  # --- Add Genome Label ---
  annoGenomeLabel(
    plot = miqtl_plot, # Attach to one of the plots, e.g., gene_plot
    params = params_genome,
    x = PLOT_PARAMS$x, y = PLOT_PARAMS$plot_y + PLOT_PARAMS$plot_height, # Below last track
    scale = "Mb", fontsize = 8,
    just = c("center", "top"), default.units = "inches")
  
  # --- Add Title ---
  plotText(
    label = paste("QTL:", locus_info$trait, "(",locus_info$drug, ") - Chr",
                  current_chr, "Peak:", round(locus_info$peak_pos, 2), "Mb"),
    x = PLOT_PARAMS$x, y = 0.1, just = c("center", "top"),
    fontface = "bold", fontsize = 9, default.units = "inches")
  
  
  dev.off()  










