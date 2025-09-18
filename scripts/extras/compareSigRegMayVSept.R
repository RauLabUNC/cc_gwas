library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(igraph, lib.loc = "/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib")

may <- read.csv("results/joinLoci/geneTables/all_significant_regions_summary_may.csv")
sept <- read.csv("data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv")

may$Peak_SNP_ID %in% sept$Peak_SNP_ID
sept$Peak_SNP_ID %in% may$Peak_SNP_ID
both <- intersect(may$Peak_SNP_ID, sept$Peak_SNP_ID)
length(both)

# Check the merged tables
trait_may <- read.csv("traitLoci.csv")
trait_sept <- read.csv("data/processed/joinLoci/relational_tables/traitLoci.csv")

str(trait_may)
str(trait_sept)

# They dont match. May has 60k rows, sept has 64 total. Sept is merged, may is not?
trait_may$Method %>% table()
trait_sept$Method %>% table()


str(may)
str(sept)
length(trait_may$Dependent_Variable) # 60k
length(trait_sept$Dependent_Variable) # 64
length(may$trait)
length(sept$trait)


### Group Loci by chromsome 
sig_regions <- sept |>
  mutate(start = pmin(upper_pos_lod_drop, lower_pos_lod_drop),
         end   = pmax(upper_pos_lod_drop, lower_pos_lod_drop))

gr <- GRanges(seqnames = paste0("chr", sig_regions$chr),        # add “chr” prefix so all seqlevels look alike
              ranges   = IRanges(sig_regions$start, sig_regions$end),
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

# --- 5. Main Loop to Generate Packets --- ####

# Define which loci to process.
group_to_process <- unique(sig_regions$locus_cluster)