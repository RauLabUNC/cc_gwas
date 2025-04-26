# This script is meant to combine CSVs of significant loci into a single CSV

# Load lib
library(data.table)
library(tidyverse)

# List the file of either type in the sig_loci directory
loci.dir <- "data/processed/sig_loci"

treatments <- c("control", "iso")
files <- lapply(treatments, function(x){
                   temp.files <- list.files(file.path(loci.dir, x))
                   temp.files <- file.path(loci.dir, x, temp.files)}) |> unlist()

# Load and bind them together
loci.dfs <- lapply(files, read.csv) |> data.table::rbindlist()

# Save
write.csv(loci.dfs, file.path(loci.dir, "allSigLoci.csv"))


# Create the plot
ggplot(loci.dfs, aes(x = -log10(threshold), fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(title = "P-values in Control vs Iso",
       x = "-log10(p-value)",
       y = "Density") +
  theme_minimal()

loci.table <- table(loci.dfs$loci)

table(loci.table)
multi.loci <- loci.table |> 
  as.data.frame() |> 
  arrange(desc(Freq)) |> 
  filter(Freq > 5) |> 
  pull(Var1)


multi.loci.df <- loci.dfs |> 
  filter(loci %in% multi.loci)

write.csv(multi.loci.df, "data/processed/sig_loci/multiSig.csv")
