# This script is meant to check for overlaps between the SNPs..
# .. sig for any trait and those for cell abundance. All results.. 
# .. are from Christoph's GWAS runs in late-March

# Load the tables
traits_iso <- read.csv("data/processed/christophGWAS/Iso_pvals.csv")
traits_ctrl <- read.csv("data/processed/christophGWAS/Ctrl_pvals.csv")
traits_delta <- read.csv("data/processed/christophGWAS/PerChange_pvals.csv")

cells <- read.csv("data/processed/christophGWAS/RoughPercentCellTypes_pvals.csv")

# load libs
library(tidyverse)

#### Data cleaning #### 
# Make a column for the smallest p value in each row
cols_cell <- colnames(cells)[6:19]
cells$p_min <- do.call(pmin, c(cells[cols_cell], na.rm = TRUE))


# run for all trait dfs
cols <- colnames(traits_iso)[6:57]
traits_iso$p_min <- do.call(pmin, c(traits_iso[cols], na.rm = TRUE))
traits_ctrl$p_min <- do.call(pmin, c(traits_ctrl[cols], na.rm = TRUE))
traits_delta$p_min <- do.call(pmin, c(traits_delta[cols], na.rm = TRUE))

## Filter to significant rows
sig_cells <- cells |> filter(p_min < 10^-3) |> ungroup() |> select(-p_min)

sig_traits_iso <- traits_iso |> filter(p_min < 10^-4)|> select(-p_min)
sig_traits_ctrl <- traits_ctrl |> filter(p_min < 10^-4) |> select(-p_min)
sig_traits_delta <- traits_delta |> filter(p_min < 10^-4) |> select(-p_min)

# Make this into a long, merged format 
sig_traits_iso <- sig_traits_iso |> 
  select(-MAF) |> 
  pivot_longer(cols = cols, names_to = "trait", values_to = "p") |> 
  filter(p < 10^-4) |> 
  mutate(context = "Iso")

sig_traits_ctrl <- sig_traits_ctrl |> 
  select(-MAF) |> 
  pivot_longer(cols = cols, names_to = "trait", values_to = "p") |> 
  filter(p < 10^-4) |> 
  mutate(context = "Control")

sig_traits_delta <- sig_traits_delta |> 
  select(-MAF) |> 
  pivot_longer(cols = cols, names_to = "trait", values_to = "p") |> 
  filter(p < 10^-4) |> 
  mutate(context = "Delta")

sig_traits <- rbind(sig_traits_ctrl, sig_traits_delta, sig_traits_iso) |> arrange(p)

# Make cells into a long format as with traits
cells_long <- cells |> 
  select(-MAF) |> 
  pivot_longer(cols = cols_cell, names_to = "trait", values_to = "p") |> 
  filter(p < 0.005)

cells_long$type <- lapply(str_split(cells_long$trait, "_"), "[[", 1) |> unlist()
cells_long$context <-  lapply(str_split(cells_long$trait, "_"), "[[", 2) |> unlist()

cells_long <- cells_long |> select(-trait, -p_min) |> arrange(p)

# Subset cells long and sig_traits to overlapping Names

cell_trait_sig <- inner_join(cells_long, sig_traits, by = c("Name", "Chr", "Pos_mm10", "Pos_mm39"))

cell_trait_sig <- cell_trait_sig |> 
  group_by(Chr) |> 
  arrange(Pos_mm10)

write.csv(cell_trait_sig, "data/processed/christophGWAS/IMGScellTraitOverlaps.csv", row.names = F)
