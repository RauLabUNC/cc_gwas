
miqtl_res <- read.csv("data/processed/sig_loci/allSigLoci.csv", row.names = 1 )
cells <- read.csv("data/processed/christophGWAS/RoughPercentCellTypes_pvals.csv")

# load libs
library(tidyverse)

#### Data cleaning #### 
# Make a column for the smallest p value in each row
cols_cell <- colnames(cells)[6:19]

cells_long <- cells |> 
  select(-MAF) |> 
  pivot_longer(cols = cols_cell, names_to = "phenotype", values_to = "p.value") |> 
  filter(p.value < 0.005) |> 
  mutate(loci = Name,
         position = Pos_mm39,
         chromosome = Chr) |> 
  select(-Name, -Chr, -Pos_mm39)

cells_long$type <- lapply(str_split(cells_long$phenotype, "_"), "[[", 1) |> unlist()
cells_long$treatment <-  lapply(str_split(cells_long$phenotype, "_"), "[[", 2) |> unlist()

cells_long <- cells_long |> select(-phenotype, -p_min) |> arrange(p.value)

# Subset cells long and sig_traits to overlapping Names

cell_trait_sig <- inner_join(cells_long, miqtl_res, by = c("loci"))
