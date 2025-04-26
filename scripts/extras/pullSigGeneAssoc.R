gene_corr <- read.csv("data/processed/christophGWAS/Gene_Trait_Cors_Controls_t21425_p31825.csv", row.names = 1)

gene_corr$gene <- row.names(gene_corr)
# Convert from wide to long format
gene_corr_long <- gene_corr %>%
  # Create trait name pattern for .R and .P suffixes
  pivot_longer(
    cols = -gene,  # Keep gene as ID column
    names_to = c("trait", ".value"),
    names_pattern = "(.*)\\.([RP])$",
    values_drop_na = TRUE
  )

sig_gene <- gene_corr_long |> 
              filter(P < 0.05)
# Check the first few rows
head(gene_corr_long)
