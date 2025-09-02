#!/usr/bin/env Rscript
# Test script for NRVM expression data analysis
# Analyzes CPM values for Itgb1, Ttn, and Myh6 genes

# NOTE: Request interactive session before running:
# srun --pty --mem=8G --time=1:00:00 --cpus-per-task=4 /bin/bash
# Then: source /proj/raulab/users/brian/claude-test/activate_miqtl.sh

# Load required libraries
# Using individual packages instead of tidyverse bundle
library(dplyr)
library(tidyr) 
library(data.table)
library(ggplot2)
library(patchwork)

# Set theme for plots
theme_set(theme_bw())

# Define paths (from the reference script)
base_path   <- "/proj/raulab/projects/cc_gwas"
nrvm_counts_path <- file.path(base_path, "data/processed/joinLoci/nrvms/bulk_gene.csv")
nrvm_meta_path <- file.path(base_path, "data/processed/joinLoci/nrvms/phenotypes.csv")

# Target genes of interest
target_genes <- c("Itgb1", "Ttn", "Myh6")

cat("====================================\n")
cat("NRVM Expression Analysis\n")
cat("====================================\n\n")

# === 1. Read data ===
cat("Reading NRVM count data...\n")
NRVM_counts <- fread(nrvm_counts_path)
setnames(NRVM_counts, 1, "mouse_gene_symbol")   # first column is gene symbol

cat("Reading metadata...\n")
NRVM_meta <- fread(nrvm_meta_path)

# Get control sample IDs
ctrl_ids <- NRVM_meta[treatment == "Ctl", sprintf("BCJ_%s", sample.id)]
cat(sprintf("Found %d control samples\n\n", length(ctrl_ids)))

# === 2. CPM transformation ===
cat("Converting to CPM...\n")
counts_mat <- as.matrix(NRVM_counts[, ..ctrl_ids])
rownames(counts_mat) <- NRVM_counts$mouse_gene_symbol

# Calculate library sizes and CPM
lib_sizes <- colSums(counts_mat)
CPM <- t(t(counts_mat) / (lib_sizes / 1e6))

# Calculate average CPM across samples
avgenes_cpm <- rowMeans(CPM)

# Create expression data table
NRVM_expr <- data.table(
  mouse_gene_symbol = NRVM_counts$mouse_gene_symbol,
  avgenes_cpm = avgenes_cpm
)

# === Create a separate plot for expression by treatment ===
cat("\nGenerating treatment-specific plot...\n")

# 1. Get all sample IDs and their treatments
all_sample_info <- NRVM_meta[, .(sample_col = sprintf("BCJ_%s", sample.id), treatment)]

# 2. Recalculate CPM for ALL samples
all_sample_ids <- all_sample_info$sample_col
counts_mat_all <- as.matrix(NRVM_counts[, ..all_sample_ids])
rownames(counts_mat_all) <- NRVM_counts$mouse_gene_symbol
lib_sizes_all <- colSums(counts_mat_all)
CPM_all <- t(t(counts_mat_all) / (lib_sizes_all / 1e6))

# 3. Prepare data for plotting
target_cpm_all <- CPM_all[rownames(CPM_all) %in% target_genes, , drop = FALSE]
if (nrow(target_cpm_all) > 0) {
  plot_data_treatment <- as.data.table(target_cpm_all, keep.rownames = "gene")
  plot_data_treatment <- melt(plot_data_treatment, 
                              id.vars = "gene", 
                              variable.name = "sample_col", 
                              value.name = "cpm")
  
  # Add treatment info
  plot_data_treatment[all_sample_info, on = "sample_col", treatment := i.treatment]
  
  # Calculate log2(CPM+1)
  plot_data_treatment[, log2_cpm := log2(cpm + 1)]
  
  # 4. Create the plot
  p_treatment <- ggplot(plot_data_treatment, aes(x = gene, y = log2_cpm, fill = treatment)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
               alpha = 0.6, size = 1.5) +
    labs(
      title = "Expression of Target Genes by Treatment",
      subtitle = sprintf("All samples (n=%d)", length(all_sample_ids)),
      x = "Gene",
      y = "log2(CPM + 1)",
      fill = "Treatment"
    ) +
    theme(legend.position = "bottom")
  
  # 5. Save the plot
  ggsave("nrvm_expression_by_treatment.pdf", p_treatment, width = 8, height = 6)
  cat("Plot saved as 'nrvm_expression_by_treatment.pdf'\n")
} else {
  cat("No target genes found for treatment-specific plot.\n")
}

# === 3. Extract target genes ===
target_data <- NRVM_expr[mouse_gene_symbol %in% target_genes]
setkey(target_data, mouse_gene_symbol)

for (gene in target_genes) {
  gene_cpm <- target_data[mouse_gene_symbol == gene, avgenes_cpm]
  if (length(gene_cpm) > 0) {
    cat(sprintf("%-10s: %12.2f CPM\n", gene, gene_cpm))
    
    # Also show individual sample values
    gene_idx <- which(rownames(CPM) == gene)
    if (length(gene_idx) > 0) {
      sample_cpms <- CPM[gene_idx, ]
      cat(sprintf("  Range: %.2f - %.2f CPM\n", min(sample_cpms), max(sample_cpms)))
      cat(sprintf("  SD: %.2f CPM\n", sd(sample_cpms)))
    }
  } else {
    cat(sprintf("%-10s: NOT FOUND\n", gene))
  }
}

# === 4. Create visualizations ===
cat("\nGenerating plots...\n")

# Prepare data for plotting
plot_data <- data.table(
  gene = NRVM_expr$mouse_gene_symbol,
  cpm = NRVM_expr$avgenes_cpm,
  log2_cpm = log2(NRVM_expr$avgenes_cpm + 1),
  is_target = NRVM_expr$mouse_gene_symbol %in% target_genes
)

# Add target gene labels
plot_data[is_target == TRUE, label := gene]

# Plot 1: Distribution of CPM values (log2 scale)
p1 <- ggplot(plot_data, aes(x = log2_cpm)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(data = plot_data[is_target == TRUE], 
             aes(xintercept = log2_cpm, color = gene),
             linetype = "dashed", size = 1) +
  geom_text(data = plot_data[is_target == TRUE],
            aes(x = log2_cpm, y = Inf, label = gene, color = gene),
            vjust = 2, hjust = -0.1, angle = 90, size = 3) +
  labs(title = "Distribution of Gene Expression (All Genes)",
       x = "log2(CPM + 1)",
       y = "Number of Genes",
       color = "Target Gene") +
  theme(legend.position = "bottom")

# Plot 2: Density plot with target genes marked
p2 <- ggplot(plot_data, aes(x = log2_cpm)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_rug(data = plot_data[is_target == TRUE], 
           aes(color = gene), sides = "b", size = 2) +
  geom_point(data = plot_data[is_target == TRUE],
             aes(y = 0, color = gene), size = 4, shape = 17) +
  labs(title = "Expression Density with Target Genes",
       x = "log2(CPM + 1)",
       y = "Density",
       color = "Target Gene") +
  theme(legend.position = "bottom")

# Plot 3: Expression rank plot
plot_data[, rank := rank(-cpm)]
plot_data[is_target == TRUE, rank_label := sprintf("%s\n(Rank: %d)", gene, rank)]

p3 <- ggplot(plot_data, aes(x = rank, y = log2_cpm)) +
  geom_line(color = "gray50", alpha = 0.5) +
  geom_point(data = plot_data[is_target == TRUE],
             aes(color = gene), size = 4) +
  geom_text(data = plot_data[is_target == TRUE],
            aes(label = rank_label, color = gene),
            vjust = -1, size = 3) +
  labs(title = "Gene Expression Rank Plot",
       x = "Gene Rank (by CPM)",
       y = "log2(CPM + 1)",
       color = "Target Gene") +
  scale_x_continuous(labels = scales::comma) +
  theme(legend.position = "bottom")

# Plot 4: Box plot of sample-level expression for target genes
# Prepare sample-level data for target genes
sample_data_list <- list()
for (gene in target_genes) {
  gene_idx <- which(rownames(CPM) == gene)
  if (length(gene_idx) > 0) {
    sample_data_list[[gene]] <- data.table(
      gene = gene,
      sample = colnames(CPM),
      cpm = CPM[gene_idx, ],
      log2_cpm = log2(CPM[gene_idx, ] + 1)
    )
  }
}

if (length(sample_data_list) > 0) {
  sample_data <- rbindlist(sample_data_list)
  
  p4 <- ggplot(sample_data, aes(x = gene, y = log2_cpm, fill = gene)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_point(position = position_jitter(width = 0.2), 
               alpha = 0.5, size = 2) +
    labs(title = "Sample-level Expression of Target Genes",
         x = "Gene",
         y = "log2(CPM + 1)",
         fill = "Gene") +
    theme(legend.position = "none")
} else {
  p4 <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "No target genes found") +
    theme_void()
}

# Combine plots
combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "NRVM Expression Analysis: Itgb1, Ttn, and Myh6",
    subtitle = sprintf("Control samples only (n=%d)", length(ctrl_ids))
  )

# Save plot
ggsave("nrvm_expression_analysis.pdf", 
       combined_plot, 
       width = 14, height = 10)

cat("\nPlot saved as 'nrvm_expression_analysis.pdf'\n")

# === 5. Summary statistics ===
cat("\n====================================\n")
cat("Summary Statistics\n")
cat("====================================\n")

# Overall expression statistics
cat("\nOverall expression distribution:\n")
cat(sprintf("  Total genes: %d\n", nrow(NRVM_expr)))
cat(sprintf("  Genes with CPM > 1: %d\n", sum(NRVM_expr$avgenes_cpm > 1)))
cat(sprintf("  Genes with CPM > 10: %d\n", sum(NRVM_expr$avgenes_cpm > 10)))
cat(sprintf("  Genes with CPM > 100: %d\n", sum(NRVM_expr$avgenes_cpm > 100)))

# Percentile information for target genes
cat("\nTarget genes percentile ranks:\n")
for (gene in target_genes) {
  gene_cpm <- target_data[mouse_gene_symbol == gene, avgenes_cpm]
  if (length(gene_cpm) > 0) {
    percentile <- ecdf(NRVM_expr$avgenes_cpm)(gene_cpm) * 100
    cat(sprintf("  %s: %.1f percentile\n", gene, percentile))
  }
}

# === 6. Export results ===
cat("\nExporting results...\n")

# Export target gene data
fwrite(target_data, "target_genes_cpm.csv")

# Export full CPM matrix for target genes (if needed for downstream analysis)
target_cpm_matrix <- CPM[rownames(CPM) %in% target_genes, , drop = FALSE]
if (nrow(target_cpm_matrix) > 0) {
  target_cpm_dt <- data.table(
    gene = rownames(target_cpm_matrix),
    as.data.table(target_cpm_matrix)
  )
  fwrite(target_cpm_dt, "target_genes_all_samples_cpm.csv")
  cat("  Exported 'target_genes_all_samples_cpm.csv'\n")
}

cat("\n====================================\n")
cat("Analysis complete!\n")
cat("====================================\n")