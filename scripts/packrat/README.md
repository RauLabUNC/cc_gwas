# PackRat: Packet-based Rating for GWAS Candidate Gene Prioritization

**Authors:** Brian Gural, Anh Luu, Todd Kimball, Christoph D. Rau
**Reference:** Gural et al. (in preparation)

## Overview

PackRat is a modular framework for creating standardized "gene packets" - evidence summaries that facilitate candidate gene prioritization from broad GWAS/QTL intervals. Each packet consolidates study-specific data (expression, eQTL, variants) with public database annotations (phenotypes, disease associations) into uniform outputs for team-based curation.

## Design Philosophy

PackRat follows a **relational data model** (see Figure B in manuscript):
- Genes, variants, phenotypes, and QTL data exist in separate tables
- Functions join these tables based on user needs
- Outputs are standardized but flexible to accommodate different study designs

## Module Organization

```
packrat/
├── join_gene_tables.R      # Join relational tables (genes + annotations)
├── phenotype_extract.R     # Customizable phenotype/disease keyword filtering
├── excel_generation.R      # Multi-sheet Excel workbook creation
├── summary_reports.R       # README and text summary generation
└── plot_locus.R            # LocusZoom-style visualization (plotgardener)
```

### 1. `join_gene_tables.R`
**Purpose:** Build master gene tables from relational data sources

**Key Functions:**
- `build_gene_table()` - Join genes with orthology, variants, expression, eQTL, phenotypes
- `get_genes_in_region()` - Extract genes within genomic coordinates
- `add_locus_membership()` - Add binary columns for multi-locus packets
- `summarize_variants()` - Collapse variant-level data to gene-level

**Flexibility:** Accepts any combination of optional data sources. Joins by gene ID or symbol.

### 2. `phenotype_extract.R`
**Purpose:** Filter genes by phenotype/disease associations using keyword matching

**Key Functions:**
- `extract_phenotypes()` - General keyword-based phenotype extraction
- `extract_human_traits()` - Human disease associations (wrapper)
- `extract_mouse_phenotypes()` - Mouse phenotype associations (wrapper)
- `format_human_only()`, `format_mouse_only()`, `format_both()` - Formatted output for reports

**Predefined Keywords:**
- `cardiac_keywords()` - Cardiovascular phenotypes
- `metabolic_keywords()` - Metabolic phenotypes
- `neurological_keywords()` - Neurological phenotypes

**Custom Keywords:** Users can define their own for any phenotype class.

### 3. `excel_generation.R`
**Purpose:** Create multi-sheet Excel workbooks for manual curation

**Key Functions:**
- `create_gene_workbook()` - Main function, highly generalizable
- `create_gene_workbook_consolidated()` - Version with consolidated phenotype sheets
- `prepare_mouse_phenotype_table()` - Format MouseMine data
- `prepare_human_disease_table()` - Format Open Targets data
- `create_variant_table()` - Export variant CSV

**Features:**
- Conditional formatting (highlight Yes/No columns)
- Auto-filtering and frozen panes
- Text wrapping and auto-width columns
- Customizable column ordering

### 4. `summary_reports.R`
**Purpose:** Generate README files for quick locus assessment

**Key Functions:**
- `create_locus_readme()` - General README with locus stats
- `create_phenotype_readme()` - README highlighting phenotype matches
- `create_top_candidates_csv()` - Export filtered candidate gene table
- `create_locus_id()` - Generate file-safe locus identifiers

**Output:** Markdown-formatted text files with gene counts, filtering criteria, and top candidates.

### 5. `plot_locus.R`
**Purpose:** LocusZoom-style visualizations (requires plotgardener)

**Key Functions:**
- `plot_locus_zoom()` - Main plotting function
- Panel-specific helpers (LOD scores, founder effects, gene tracks)

**Current Implementation:** Mouse-specific (mm39), uses plotgardener and TxDb/OrgDb packages.

## Example Usage

### Minimal Example: Gene List → Excel Packet

```r
# Load modules
source("scripts/packrat/join_gene_tables.R")
source("scripts/packrat/excel_generation.R")

# Define genes in your locus
genes <- data.frame(
  gene_symbol = c("Abcb10", "Acsl5", "Pdlim5"),
  ensembl_id = c("ENSMUSG001", "ENSMUSG002", "ENSMUSG003"),
  chr = c(8, 8, 8),
  start_bp = c(28300000, 30500000, 35200000),
  end_bp = c(28350000, 30550000, 35250000)
)

# Create Excel workbook
create_gene_workbook(
  gene_table = genes,
  output_file = "my_locus_genes.xlsx",
  essential_cols = c("gene_symbol", "ensembl_id", "chr")
)
```

### With Expression and eQTL Data

```r
# Add expression data
expression <- data.frame(
  gene_symbol = c("Abcb10", "Acsl5", "Pdlim5"),
  avg_cpm_ctrl = c(150, 45, 220),
  log2fc_treatment = c(0.5, -0.2, 1.2),
  padj = c(0.001, 0.8, 0.0001)
)

# Add eQTL status
eqtl <- data.frame(
  gene_symbol = c("Abcb10", "Pdlim5"),
  has_cis_eqtl = "Yes",
  eqtl_lod = c(8.5, 12.3)
)

# Build master table
master_table <- build_gene_table(
  genes_in_locus = genes,
  gene_symbol_col = "gene_symbol",
  expression = expression,
  eqtl = eqtl
)

# Create workbook
create_gene_workbook(
  gene_table = master_table,
  output_file = "my_locus_with_expression.xlsx",
  highlight_cols = c("has_cis_eqtl")
)
```

### With Phenotype Filtering

```r
source("scripts/packrat/phenotype_extract.R")

# Add phenotype data (from MouseMine or Open Targets)
phenotypes <- data.frame(
  gene_symbol = c("Abcb10", "Abcb10", "Pdlim5"),
  phenotype = c("abnormal heart morphology",
                "increased response to stress",
                "dilated cardiomyopathy")
)

# Filter for cardiac genes
cardiac_keywords <- cardiac_keywords()  # Use predefined keywords
cardiac_genes <- genes %>%
  filter(gene_symbol %in% phenotypes$gene_symbol)

# Create summary
create_phenotype_readme(
  locus_info = list(chr = 8, start = 28000000, end = 41000000),
  gene_table = master_table,
  phenotype_keywords = cardiac_keywords,
  output_file = "README_cardiac_genes.txt"
)
```

## Adapting for Your Study

### 1. Different Organisms
- Modify `join_gene_tables.R`: change default column names (e.g., `human_ensembl_id`)
- Modify `plot_locus.R`: use appropriate genome assembly and TxDb/OrgDb

### 2. Different Phenotype Databases
- Create custom keyword sets in `phenotype_extract.R`
- Use `extract_phenotypes()` with your delimiter and column names

### 3. Additional Data Sources
- Use `custom_data` parameter in `build_gene_table()`
- Tables just need a matching gene ID or symbol column

### 4. Custom Excel Layouts
- Use `essential_cols` to reorder columns
- Use `highlight_cols` to add conditional formatting
- Create custom preparation functions for your data format

## Integration with Main Pipeline

The refactored `21_makeLociPackets.R` script loads these modules and orchestrates:
1. Loading study-specific data (QTL scans, expression, variants)
2. Clustering overlapping loci
3. For each cluster:
   - Get genes in region
   - Build master table (join annotations)
   - Generate Excel workbook
   - Generate plots
   - Generate README

See `scripts/packets/21_makeLociPackets.R` for the full implementation.

## Dependencies

**Required:**
- `data.table` - Fast data manipulation
- `dplyr` - Data wrangling
- `openxlsx` - Excel file creation

**Optional (for plotting):**
- `plotgardener` - LocusZoom-style plots
- `GenomicRanges` - Genomic interval operations
- `TxDb.*` and `org.*.eg.db` - Organism annotations

## Citation

If you use PackRat in your research, please cite:

Gural et al. (in preparation). "______: a Semi-Automated Framework for Prioritizing Candidate Genes from Broad GWAS Intervals." *BMC Research Notes*.

## Contact

For questions or contributions:
- Brian Gural: [email]
- Christoph Rau: [email]
- GitHub: https://github.com/RauLabUNC/cc_gwas
