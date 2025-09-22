# Load
data <- read.csv("all_significant_regions_summary.csv", stringsAsFactors = FALSE)

# Y: traits with ≥1 significant locus; Z: traits with hits in only one treatment
Y <- length(unique(data$trait))
Z <- sum(tapply(data$drug, data$trait, function(x) length(unique(x)) == 1))

# n: total loci
n <- nrow(data)

# Overlap logic (by chr, using interval intersection)
start <- pmin(data$lower_pos_lod_drop, data$upper_pos_lod_drop)
end   <- pmax(data$lower_pos_lod_drop, data$upper_pos_lod_drop)
chr_eq <- outer(data$chr, data$chr, "==")
ov <- outer(start, end, function(s1, e2) s1 <= e2) & outer(end, start, function(e1, s2) e1 >= s2)
ov_mat <- chr_eq & ov; diag(ov_mat) <- FALSE

# N: loci that overlap at least one other locus
N_overlap <- sum(rowSums(ov_mat) > 0)

# Consolidate to multi-trait loci and count chromosomes (no igraph)
# Build edge list from overlap matrix (upper triangle to avoid duplicates)
edges <- which(ov_mat & upper.tri(ov_mat), arr.ind = TRUE)

# Union–Find (Disjoint Set Union) to get connected components
parent <- seq_len(nrow(data))
find <- function(x) {
  while (parent[x] != x) {
    parent[x] <<- parent[parent[x]]  # path compression
    x <- parent[x]
  }
  x
}
unite <- function(a, b) {
  ra <- find(a); rb <- find(b)
  if (ra != rb) parent[rb] <<- ra
}
if (nrow(edges) > 0) {
  apply(edges, 1, function(e) unite(e[1], e[2]))
}
comp <- vapply(seq_along(parent), find, integer(1))

library(dplyr)
by_comp <- data %>%
  mutate(component = comp) %>%
  group_by(component)

traits_per_cluster <- by_comp %>%
  summarise(n_traits = n_distinct(trait), .groups = "drop")

multi_components <- traits_per_cluster %>%
  filter(n_traits >= 2) %>%
  pull(component)

N_multi_trait_loci <- length(multi_components)

chroms_per_cluster <- by_comp %>%
  filter(component %in% multi_components) %>%
  summarise(chrs = list(unique(chr)), .groups = "drop") %>%
  pull(chrs)

N_chromosomes <- length(unique(unlist(chroms_per_cluster)))

# Optional: view all in one vector
c(Y = Y, Z = Z, n = n, N_overlap = N_overlap,
  N_multi_trait_loci = N_multi_trait_loci, N_chromosomes = N_chromosomes)
