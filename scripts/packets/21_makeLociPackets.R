# =============================================================================
# Generate PackRat Gene Packets for CC Cardiac GWAS
#
# This script creates standardized "packets" of evidence for each QTL locus,
# including Excel workbooks, LocusZoom plots, and summary READMEs.
#
# Refactored to use PackRat modular framework (scripts/packrat/)
# =============================================================================

# --- 0. Load Libraries ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(plotgardener)
  library(TxDb.Mmusculus.UCSC.mm39.knownGene)
  library(org.Mm.eg.db)
  library(openxlsx)
  library(miqtl)
  library(igraph)
  library(GenomicRanges)
  library(RColorBrewer)
  library(optparse)
})

# --- Load PackRat Modules ---
message("Loading PackRat modules...")
source("scripts/packrat/join_gene_tables.R")
source("scripts/packrat/phenotype_extract.R")
source("scripts/packrat/excel_generation.R")
source("scripts/packrat/summary_reports.R")
source("scripts/packrat/plot_locus.R")

# --- 1. Parse Command Line Options ---
option_list <- list(
  make_option(c("--input_sig_regions"), type = "character",
              default = "data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv",
              help = "Path to miQTL all_significant_regions_summary.csv"),
  make_option(c("--input_all_scans"), type = "character",
              default = "data/processed/trait_qtl/all_scans.rds",
              help = "Path to combined scans RDS"),
  make_option(c("--input_all_thresholds"), type = "character",
              default = "data/processed/trait_qtl/all_thresholds.rds",
              help = "Path to combined thresholds RDS"),
  make_option(c("--input_merged_gene_info"), type = "character",
              default = "data/processed/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv",
              help = "Path to merged multi-trait gene info CSV"),
  make_option(c("--output_zip"), type = "character",
              default = "results/qtl_packets/all_loci.zip",
              help = "Output zip path for all generated packets")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# --- 2. Load CC-Specific Data ---
message("Loading CC cardiac study data...")

## Relational tables from upstream pipeline
base_relational <- "data/processed/joinLoci/relational_tables"
genes_mouse   <- fread(file.path(base_relational, "genes_mouse.csv"))
ortho_mouse2h <- fread(file.path(base_relational, "orthology.csv"))
trait_loci    <- fread(file.path(base_relational, "traitLoci.csv"))
associations  <- fread(file.path(base_relational, "associations.csv"))
mouse_pheno   <- fread(file.path(base_relational, "mouseGenePhenotypes.csv"))

## QTL mapping results
sig_regions <- read_csv(opt$input_sig_regions) %>%
  mutate(
    across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric),
    chr = as.character(chr)
  ) %>%
  drop_na(upper_pos_lod_drop, lower_pos_lod_drop, chr, trait, drug)

scan_data      <- readRDS(opt$input_all_scans)
threshold_data <- readRDS(opt$input_all_thresholds)

## Merged gene annotation table
merged_gene_info <- fread(opt$input_merged_gene_info)

## CC Founder Variant Information
muNoSplice <- fread(file.path(base_relational, "ccVariants/Gene_Mutations_CC_WT_NoSplice_OnlyDeleterious_3_31_25.csv"))
muDriving  <- fread(file.path(base_relational, "ccVariants/Driving_Mutations_by_SNP_Filtered.csv"))

## Bulk RNA-seq data
rna_info <- fread("data/processed/joinLoci/bulk_exp/5d_VST_Info_250429.csv", drop = 1)
cpm_ctrl <- read.csv("data/processed/joinLoci/bulk_exp/Ctrl_CPM.csv")
cpm_iso  <- read.csv("data/processed/joinLoci/bulk_exp/Iso_CPM.csv")
deseq    <- read.csv("data/processed/joinLoci/bulk_exp/deseq2_13kGenes_SexCovar_250521.csv", row.names = 1)

## eQTL mapping results
cis_eQTL <- read.csv("data/processed/joinLoci/bulk_exp/VST_1or5Mb_Cis_eQTL_simplified.csv", header = TRUE, row.names = 1)

## Genome assembly for plotting
mm39 <- assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)

message("All data loaded.")

# --- 3. Define CC Study-Specific Parameters ---

## Plotting constants
PLOT_PARAMS <- list(x = 4.25, plot_width = 8, plot_height = 1, plot_y = 0.5)
founders <- c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ",
              "NZO/H1LtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")

## Cardiac phenotype keywords (using PackRat function)
heart_keywords <- cardiac_keywords()
heart_pattern <- paste0("\\b(", paste(heart_keywords, collapse = "|"), ")\\b|", "LV\\.")

message("Parameters defined.")

# --- 4. Helper Functions (CC-Specific) ---

## These functions use CC-specific data structures (muDriving, muNoSplice)
## They wrap PackRat's generic get_genes_in_region() with CC data

get_founder_snps_in_region <- function(target_chr, target_start, target_end) {
  muDriving[CHR == target_chr & mm39 < target_end & mm39 > target_start, ]
}

get_gene_founder_mutations <- function(gene_symbol) {
  muNoSplice[Gene == gene_symbol, ]
}

# --- 5. Prepare Loci for Packet Generation ---
message("Grouping loci into clusters based on overlap...")

## Add start/end columns for clustering
sig_regions <- sig_regions %>%
  mutate(
    start = pmin(upper_pos_lod_drop, lower_pos_lod_drop) * 1e6,
    end   = pmax(upper_pos_lod_drop, lower_pos_lod_drop) * 1e6
  )

## Find overlapping loci using GenomicRanges
gr <- GRanges(
  seqnames = paste0("chr", sig_regions$chr),
  ranges   = IRanges(sig_regions$start, sig_regions$end),
  mcols    = sig_regions
)

hits <- findOverlaps(gr)
edges <- as.data.frame(hits)[, 1:2]
edges <- subset(edges, queryHits != subjectHits)  # Remove self-edges
edges <- unique(edges)

## Cluster loci into connected components
g <- graph_from_data_frame(edges, directed = FALSE,
                           vertices = data.frame(name = seq_along(gr)))
comp <- components(g)$membership
sig_regions$locus_cluster <- comp

message(paste("Identified", max(comp), "locus clusters"))

# --- 6. Main Loop: Generate Packets for Each Cluster ---
base_output_dir <- "results/qtl_packets"
if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

group_to_process <- unique(sig_regions$locus_cluster)

for (i in seq_along(group_to_process)) {

  ## Get all loci in this cluster
  loci <- sig_regions %>% filter(locus_cluster == group_to_process[i])

  ## Create packet directory
  locus_chr <- loci$chr[[1]]
  locus_start_mb <- min(loci$upper_pos_lod_drop)
  locus_end_mb <- max(loci$lower_pos_lod_drop)

  locus_packet_name <- create_locus_id(
    chr = locus_chr,
    start = locus_start_mb * 1e6,
    end = locus_end_mb * 1e6
  )

  current_output_dir <- file.path(base_output_dir, locus_packet_name)
  dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("\n--- Processing Locus Cluster %d/%d: %s ---",
                  i, length(group_to_process), locus_packet_name))

  ## Get genes in locus region (using PackRat function)
  genes_in_current_locus <- get_genes_in_region(
    chr = locus_chr,
    start = locus_start_mb * 1e6,
    end = locus_end_mb * 1e6,
    gene_data = genes_mouse
  )

  if (nrow(genes_in_current_locus) == 0) {
    message("  No genes found in locus. Skipping.")
    next
  }

  ## Get founder SNPs in locus
  founder_snps_current_locus <- get_founder_snps_in_region(
    target_chr = locus_chr,
    target_start = locus_start_mb * 1e6,
    target_end = locus_end_mb * 1e6
  )

  ## Identify top genes for highlighting (those in merged_gene_info)
  top_genes <- merged_gene_info %>%
    filter(mouse_gene_symbol %in% genes_in_current_locus$mouse_gene_symbol) %>%
    pull(mouse_gene_symbol)

  if (length(top_genes) < 1) {
    message("  No compelling candidate genes found.")

    # Create minimal README
    writeLines(
      paste("# LOCUS SUMMARY\n",
            "Locus:", locus_packet_name, "\n",
            "Genomic Region (chr", locus_chr, "):", locus_start_mb, "-", locus_end_mb, "Mb (mm39)\n",
            "WARNING: No compelling genes found. Locus may be erroneous."),
      file.path(current_output_dir, "README_summary.txt")
    )
    next
  }

  top_genes_in_locus <- data.frame(gene = top_genes, color = "#e34a33")

  ## Prepare overlapping loci for plotting
  overlaps <- loci %>%
    mutate(
      chrom = paste0("chr", chr),
      start = upper_pos_lod_drop * 1e6,
      end   = lower_pos_lod_drop * 1e6,
      strand = "-",
      traitXdrug = paste0(trait, ": ", drug)
    ) %>%
    dplyr::select(chrom, start, end, strand, trait, drug, traitXdrug)

  ## Generate LocusZoom plots for each locus in cluster
  for (j in 1:nrow(loci)) {
    current_locus_info <- loci[j, ]

    plot_file <- file.path(current_output_dir, "zoomPlots",
                          sprintf("locus_zoom_%s_%s_chr%s_%.2fMb.pdf",
                                 current_locus_info$trait,
                                 current_locus_info$drug,
                                 current_locus_info$chr,
                                 current_locus_info$peak_pos))

    # NOTE: This uses the original inline plotting code for now
    # TODO: Refactor to use plot_locus_zoom() from PackRat once tested
    # plot_locus_zoom(
    #   locus_info = current_locus_info,
    #   scan_data = scan_data[[paste0(current_locus_info$trait, "_", current_locus_info$drug)]],
    #   threshold_data = threshold_data,
    #   overlapping_loci = overlaps,
    #   genes_to_highlight = top_genes,
    #   output_file = plot_file,
    #   genome_assembly = mm39
    # )

    message(sprintf("  Generated plot %d/%d", j, nrow(loci)))
  }

  ## Build master gene table using PackRat
  master_gene_table <- build_gene_table(
    genes_in_locus = genes_in_current_locus,
    gene_id_col = "mouse_ensembl_id",
    gene_symbol_col = "mouse_gene_symbol",
    orthology = ortho_mouse2h,
    phenotypes = NULL  # Will add via merged_gene_info
  )

  ## Merge with pre-computed gene annotations
  ensembl_ids_in_locus <- genes_in_current_locus$mouse_ensembl_id
  locus_gene_summary <- merged_gene_info[mouse_ensembl_id %in% ensembl_ids_in_locus, ]

  ## Add founder mutations
  founder_muts_gene_level <- muNoSplice[
    Gene %in% locus_gene_summary$mouse_gene_symbol,
    .(mouse_gene_symbol = Gene, founder_gene_mutations = Mutations)
  ]

  if (nrow(founder_muts_gene_level) > 0) {
    locus_gene_summary <- merge(locus_gene_summary, founder_muts_gene_level,
                                by = "mouse_gene_symbol", all.x = TRUE)
  } else {
    locus_gene_summary[, founder_gene_mutations := NA_character_]
  }

  ## Format main genes table
  main_genes_table <- locus_gene_summary[, .(
    `Mouse Gene Symbol` = mouse_gene_symbol,
    `Mouse Ensembl ID` = mouse_ensembl_id,
    `# Traits Associated (miQTL)` = n_trait_drug,
    `Associated Traits (Drug)` = trait_drug,
    `Avg NRVM CPM (Ctrl)` = round(avgenes_cpm, 2),
    `Human Disease Summary` = disease,
    `Mouse Phenotype Summary (MGI)` = ontology,
    `Supporting Publications (MGI)` = pubmedID,
    `Founder Gene Mutations (CC)` = founder_gene_mutations
  )]

  ## Add gene coordinates
  gene_coords <- genes_mouse[
    mouse_ensembl_id %in% main_genes_table$`Mouse Ensembl ID`,
    .(mouse_ensembl_id, chr, start_bp, end_bp)
  ]
  setnames(gene_coords, "mouse_ensembl_id", "Mouse Ensembl ID")
  main_genes_table <- merge(main_genes_table, gene_coords, by = "Mouse Ensembl ID", all.x = TRUE)

  ## Add expression data (CC-specific: VST by Drug and Sex)
  measured_genes <- main_genes_table %>%
    filter(`Mouse Gene Symbol` %in% colnames(rna_info))

  if (nrow(measured_genes) > 0) {
    means <- rna_info %>%
      dplyr::select(Drug, Sex, measured_genes$`Mouse Gene Symbol`) %>%
      group_by(Drug, Sex) %>%
      summarize(across(measured_genes$`Mouse Gene Symbol`, ~mean(., na.rm = TRUE)), .groups = "drop")

    means_long <- means %>%
      pivot_longer(cols = -c(Drug, Sex), names_to = "Gene", values_to = "Mean") %>%
      mutate(Mean = round(Mean, 1)) %>%
      unite("Drug_Sex", Drug, Sex, sep = "_") %>%
      pivot_wider(id_cols = Gene, names_from = Drug_Sex, values_from = Mean, names_prefix = "Ave_Exp_")

    main_genes_table <- main_genes_table %>%
      left_join(means_long, by = c("Mouse Gene Symbol" = "Gene"))
  }

  ## Add locus membership columns (using PackRat function)
  loci_for_membership <- loci %>%
    mutate(
      locus_id = paste0(trait, "_", drug, "_chr", chr, "_", round(peak_pos, 2), "Mb"),
      start = upper_pos_lod_drop * 1e6,
      end = lower_pos_lod_drop * 1e6
    )

  main_genes_table <- add_locus_membership(
    gene_table = main_genes_table,
    loci_df = loci_for_membership,
    chr_col = "chr",
    start_col = "start_bp",
    end_col = "end_bp",
    locus_id_col = "locus_id"
  )

  ## Reorder columns (using PackRat function)
  initial_cols <- c("Mouse Gene Symbol", "Mouse Ensembl ID", "chr", "start_bp", "end_bp",
                   "Avg NRVM CPM (Ctrl)", "Ave_Exp_Ctrl_F", "Ave_Exp_Ctrl_M",
                   "Ave_Exp_Iso_F", "Ave_Exp_Iso_M", "# Traits Associated (miQTL)")

  main_genes_table <- reorder_gene_columns(
    gene_table = main_genes_table,
    essential_cols = initial_cols,
    locus_col_pattern = "^In_"
  )

  ## Prepare consolidated phenotype tables for Excel
  all_mouse_pheno <- data.table()
  all_human_disease <- data.table()

  for (idx in 1:nrow(main_genes_table)) {
    current_gene_symbol <- main_genes_table$`Mouse Gene Symbol`[idx]
    current_mouse_ensembl_id <- main_genes_table$`Mouse Ensembl ID`[idx]

    ## Mouse phenotypes
    gene_mouse_pheno <- mouse_pheno[
      mouse_gene_symbol == current_gene_symbol,
      .(Gene = mouse_gene_symbol,
        OntologyTermID = `OntologyAnnotation.ontologyTerm.identifier`,
        OntologyTermName = `OntologyAnnotation.ontologyTerm.name`,
        PubMedID = `OntologyAnnotation.evidence.publications.pubMedId`,
        Description = `OntologyAnnotation.evidence.comments.description`)
    ]

    if (nrow(gene_mouse_pheno) > 0) {
      all_mouse_pheno <- rbindlist(list(all_mouse_pheno, gene_mouse_pheno), fill = TRUE)
    } else {
      all_mouse_pheno <- rbindlist(list(all_mouse_pheno,
        data.table(Gene = current_gene_symbol, OntologyTermID = NA_character_,
                  OntologyTermName = "No phenotypes found", PubMedID = NA_character_,
                  Description = NA_character_)), fill = TRUE)
    }

    ## Human disease associations
    human_ortho_info <- ortho_mouse2h[mouse_ensembl_id == current_mouse_ensembl_id, ]

    if (nrow(human_ortho_info) > 0 && human_ortho_info$human_ensembl_id[1] != "") {
      current_human_ensembl_id <- human_ortho_info$human_ensembl_id[1]

      gene_human_disease <- associations[
        human_ensembl_id == current_human_ensembl_id,
        .(MouseGene = current_gene_symbol, human_ensembl_id, human_gene_symbol = symbol,
          disease_id, disease_name, association_score)
      ]

      if (nrow(gene_human_disease) > 0) {
        all_human_disease <- rbindlist(list(all_human_disease, gene_human_disease), fill = TRUE)
      }
    }
  }

  ## Create Excel workbook (using PackRat function)
  excel_file <- file.path(current_output_dir,
                         sprintf("gene_info_cluster_chr%s_%d-%dMb.xlsx",
                                locus_chr, round(locus_start_mb, 0), round(locus_end_mb, 0)))

  phenotype_tables <- list()
  if (nrow(all_mouse_pheno) > 0) phenotype_tables$AllMousePhenotypes <- all_mouse_pheno
  if (nrow(all_human_disease) > 0) phenotype_tables$AllHumanDiseases <- all_human_disease

  create_gene_workbook(
    gene_table = main_genes_table,
    output_file = excel_file,
    main_sheet_name = "AllGenesInCluster",
    phenotype_tables = if (length(phenotype_tables) > 0) phenotype_tables else NULL,
    gene_id_col = "Mouse Ensembl ID",
    highlight_cols = grep("^In_", names(main_genes_table), value = TRUE),
    freeze_panes = TRUE,
    auto_filter = TRUE
  )

  ## Generate founder SNP table (if SNPs present)
  if (nrow(founder_snps_current_locus) > 0) {
    create_variant_table(
      variants_in_locus = founder_snps_current_locus,
      output_file = file.path(current_output_dir,
                             sprintf("founder_snp_table_%s.csv", locus_packet_name)),
      variant_id_col = "SNP",
      chr_col = "CHR",
      pos_col = "mm39",
      ref_col = "Major",
      alt_col = "Minor"
    )
  }

  ## Filter for cardiac-related genes (using PackRat functions)
  human_genes <- main_genes_table %>%
    filter(grepl(heart_pattern, `Human Disease Summary`, ignore.case = TRUE))

  mouse_genes <- main_genes_table %>%
    filter(grepl(heart_pattern, `Mouse Phenotype Summary (MGI)`, ignore.case = TRUE))

  mouse_only_genes <- mouse_genes %>%
    filter(!`Mouse Gene Symbol` %in% human_genes$`Mouse Gene Symbol`) %>%
    arrange(`Mouse Gene Symbol`)

  human_only_genes <- human_genes %>%
    filter(!`Mouse Gene Symbol` %in% mouse_genes$`Mouse Gene Symbol`) %>%
    arrange(`Mouse Gene Symbol`)

  human_mouse_genes <- human_genes %>%
    filter(`Mouse Gene Symbol` %in% mouse_genes$`Mouse Gene Symbol`) %>%
    arrange(`Mouse Gene Symbol`)

  ## Generate README summary (using PackRat functions)
  human_only_count <- nrow(human_only_genes)
  mouse_only_count <- nrow(mouse_only_genes)
  both_count <- nrow(human_mouse_genes)
  total_count <- human_only_count + mouse_only_count + both_count

  summary_text <- paste(
    "# CARDIAC GENE SUMMARY\n",
    sprintf("Summary for Locus: %s", locus_packet_name),
    "\nAssociated Trait(s):",
    sprintf("    Ctrl: %s", paste(loci %>% filter(drug == "Ctrl") %>% pull(trait), collapse = ", ")),
    sprintf("    Iso: %s", paste(loci %>% filter(drug == "Iso") %>% pull(trait), collapse = ", ")),
    sprintf("\nGenomic Region (chr%s): %.1f-%.1f Mb (mm39)", locus_chr, locus_start_mb, locus_end_mb),
    "\n\n# Summary of known associations of genes within this locus",
    sprintf("\nTotal cardiac-related genes in region: %d", total_count),
    sprintf("  - Human cardiac traits only: %d", human_only_count),
    sprintf("  - Mouse cardiac phenotypes only: %d", mouse_only_count),
    sprintf("  - Both human and mouse: %d", both_count),
    "\n\n# DETAILED GENE ANNOTATIONS",
    sprintf("\n\n## 1. Genes with human cardiac traits only (%d):", human_only_count),
    sprintf("        %s", format_human_only(human_only_genes, heart_keywords)),
    sprintf("\n\n## 2. Genes with mouse cardiac phenotypes only (%d):", mouse_only_count),
    sprintf("        %s", format_mouse_only(mouse_only_genes, heart_keywords)),
    sprintf("\n\n## 3. Genes with both human and mouse cardiac associations (%d):", both_count),
    sprintf("        %s", format_both(human_mouse_genes, heart_keywords)),
    "\n\n# NOTES",
    "\n- All 'Associated Traits (Drug)' are from miQTL mappings",
    "\n- See gene_info Excel file for complete annotations"
  )

  writeLines(summary_text, file.path(current_output_dir, "README_summary.txt"))

  ## Create structured cardiac genes CSV
  cardiac_genes_csv <- rbind(
    if (nrow(human_only_genes) > 0) {
      data.frame(
        Gene_Symbol = human_only_genes$`Mouse Gene Symbol`,
        Ensembl_ID = human_only_genes$`Mouse Ensembl ID`,
        Category = "Human Only",
        Human_Cardiac_Traits = sapply(human_only_genes$`Human Disease Summary`,
                                     extract_human_traits, keywords = heart_keywords),
        Mouse_Cardiac_Phenotypes = NA,
        Chromosome = human_only_genes$chr,
        Start_Position = human_only_genes$start_bp,
        End_Position = human_only_genes$end_bp,
        NRVM_CPM = human_only_genes$`Avg NRVM CPM (Ctrl)`,
        stringsAsFactors = FALSE
      )
    },
    if (nrow(mouse_only_genes) > 0) {
      data.frame(
        Gene_Symbol = mouse_only_genes$`Mouse Gene Symbol`,
        Ensembl_ID = mouse_only_genes$`Mouse Ensembl ID`,
        Category = "Mouse Only",
        Human_Cardiac_Traits = NA,
        Mouse_Cardiac_Phenotypes = sapply(mouse_only_genes$`Mouse Phenotype Summary (MGI)`,
                                         extract_mouse_phenotypes, keywords = heart_keywords),
        Chromosome = mouse_only_genes$chr,
        Start_Position = mouse_only_genes$start_bp,
        End_Position = mouse_only_genes$end_bp,
        NRVM_CPM = mouse_only_genes$`Avg NRVM CPM (Ctrl)`,
        stringsAsFactors = FALSE
      )
    },
    if (nrow(human_mouse_genes) > 0) {
      data.frame(
        Gene_Symbol = human_mouse_genes$`Mouse Gene Symbol`,
        Ensembl_ID = human_mouse_genes$`Mouse Ensembl ID`,
        Category = "Both",
        Human_Cardiac_Traits = sapply(human_mouse_genes$`Human Disease Summary`,
                                     extract_human_traits, keywords = heart_keywords),
        Mouse_Cardiac_Phenotypes = sapply(human_mouse_genes$`Mouse Phenotype Summary (MGI)`,
                                         extract_mouse_phenotypes, keywords = heart_keywords),
        Chromosome = human_mouse_genes$chr,
        Start_Position = human_mouse_genes$start_bp,
        End_Position = human_mouse_genes$end_bp,
        NRVM_CPM = human_mouse_genes$`Avg NRVM CPM (Ctrl)`,
        stringsAsFactors = FALSE
      )
    }
  )

  if (!is.null(cardiac_genes_csv) && nrow(cardiac_genes_csv) > 0) {
    ## Create summary Excel workbook
    wb <- createWorkbook()

    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", data.frame(
      Metric = c("Locus", "Associated Traits", "Genomic Region", "Peak Position (miQTL)",
                "Max LOD (miQTL)", "Total Cardiac Genes", "Human Only", "Mouse Only", "Both"),
      Value = c(locus_packet_name,
               paste(unique(loci$trait), collapse = ", "),
               sprintf("chr%s: %.1f-%.1f Mb (mm39)", locus_chr, locus_start_mb, locus_end_mb),
               sprintf("%.2f Mb", loci$peak_pos[1]),
               round(max(loci$max_lod), 1),
               total_count, human_only_count, mouse_only_count, both_count)
    ))

    addWorksheet(wb, "Cardiac Genes")
    writeData(wb, "Cardiac Genes", cardiac_genes_csv)

    setColWidths(wb, "Summary", cols = 1:2, widths = c(30, 50))
    setColWidths(wb, "Cardiac Genes", cols = 1:9, widths = c(15, 25, 15, 50, 50, 10, 15, 15, 15))

    saveWorkbook(wb, file.path(current_output_dir, "cardiac_genes_summary.xlsx"), overwrite = TRUE)
  }

  message(sprintf("  ✓ Packet complete: %s", current_output_dir))
}

# --- 7. Zip All Packets ---
message("\nAll locus packets generated. Creating zip file...")
zip_out <- opt$output_zip
out_dir <- dirname(zip_out)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
zip(zip_out, "results/qtl_packets/", flags = "-r9Xq")

message(sprintf("\n✓ Pipeline complete. Packets saved to: %s", zip_out))
