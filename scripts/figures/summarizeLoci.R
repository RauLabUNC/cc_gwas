df <- data.frame(a=1:3, b=letters[1:3])
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(igraph)

sig_regions <- read.csv( "data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv")
genes <- read.csv("data/processed/joinLoci/relational_tables/genes_mouse.csv")
annot_genes <- read.csv("data/processed/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv")

sig_regions <- sig_regions |>
  mutate( across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
        chr = as.character(chr),
        start_bp = as.integer(floor(upper_pos_lod_drop * 1e6)),
        end_bp   = as.integer(floor(lower_pos_lod_drop * 1e6))) |>
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug) |>
  dplyr::select(-upper_pos_lod_drop, -lower_pos_lod_drop) |>
  arrange(chr, start_bp, end_bp)

seqs <- ifelse(grepl("^chr", sig_regions$chr), sig_regions$chr, paste0("chr", sig_regions$chr))
gr <- GRanges(seqnames = seqs,
              ranges   = IRanges(sig_regions$start_bp, sig_regions$end_bp),
              strand   = "*",
              mcols    = sig_regions)
hits <- GenomicRanges::findOverlaps(gr) 

edges <- as.data.frame(hits)[ , 1:2]               # queryHits | subjectHits

edges <- subset(edges, queryHits != subjectHits)   # drop the trivial self‑edges
edges <- unique(edges)                             # remove duplicated pairs

g <- graph_from_data_frame(edges,
                           directed = FALSE,
                           vertices = data.frame(name = seq_along(gr)))
## Connected components = overlapping loci 
comp <- components(g)$membership
sig_regions$locus_cluster <- comp

clean_regions <- sig_regions |> 
    mutate(size = end_bp - start_bp) |>
    arrange(size) 

summary(clean_regions$size)

# Make GRanges for genes (ensure chr prefix consistency)
gene_seqs <- ifelse(grepl("^chr", genes$chr), genes$chr, paste0("chr", genes$chr))
gr_genes <- GRanges(seqnames = gene_seqs,
                    ranges   = IRanges(genes$start_bp, genes$end_bp),
                    strand   = "*",
                    mcols    = genes)

# Overlaps (query = loci regions, subject = genes)
ov_rg <- findOverlaps(gr, gr_genes)

ov_df <- as.data.frame(ov_rg) |>
  dplyr::mutate(
    region_index = queryHits,
    gene_index   = subjectHits,
    mouse_ensembl_id = genes$mouse_ensembl_id[gene_index],
    mouse_gene_symbol = genes$mouse_gene_symbol[gene_index]
  )

# Per-region gene counts (any overlap)
region_gene_counts <- ov_df |>
    dplyr::group_by(region_index) |>
    dplyr::summarise(
        n_genes_region = dplyr::n_distinct(mouse_ensembl_id),
        genes_region_symbols = paste(sort(unique(mouse_gene_symbol)), collapse = ";"),
        genes_region_ensembl = paste(sort(unique(mouse_ensembl_id)), collapse = ";"),
        .groups = "drop"
    )

# Attach to sig_regions (fill zeros where no overlap)
sig_regions <- sig_regions |>
    dplyr::mutate(region_index = dplyr::row_number()) |>
    dplyr::left_join(region_gene_counts, by = "region_index") |>
    dplyr::mutate(
        n_genes_region = dplyr::coalesce(n_genes_region, 0L),
        genes_region_symbols = dplyr::coalesce(genes_region_symbols, ""),
        genes_region_ensembl = dplyr::coalesce(genes_region_ensembl, "")
    )

# Get unique genes across all regions
all_genes_ensembl <- sig_regions |>
    dplyr::filter(genes_region_ensembl != "") |>
    dplyr::pull(genes_region_ensembl) |>
    strsplit(";") |>
    unlist() |>
    unique()

n_unique_genes_total <- length(all_genes_ensembl)

str(annot_genes)










#write.csv(viewable, "results/sig_regions/locus_clusters.csv", row.names = FALSE)

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


# Line used to retrieve genes within a locus during packet generation
gene_data[chr == target_chr & start_bp < target_end & end_bp > target_start, ]