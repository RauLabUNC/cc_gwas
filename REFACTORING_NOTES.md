# PackRat Refactoring Notes

## What Changed

This branch refactors the CC cardiac GWAS gene packet generation pipeline to use the modular **PackRat** framework.

### New Structure

**Before:**
```
scripts/packets/21_makeLociPackets.R    # 1130 lines, monolithic
scripts/helpers/plots.R                  # Inline plotting functions
scripts/helpers/tables.R                 # Inline table generation functions
```

**After:**
```
scripts/packrat/                         # New modular framework
├── join_gene_tables.R                   # Relational table operations
├── phenotype_extract.R                  # Keyword-based filtering
├── excel_generation.R                   # Multi-sheet workbook creation
├── summary_reports.R                    # README generation
├── plot_locus.R                         # LocusZoom visualization
└── README.md                            # Usage documentation

scripts/packets/21_makeLociPackets.R     # Refactored to ~650 lines
scripts/packets/21_makeLociPackets_ORIGINAL.R  # Backup of original

scripts/helpers/
├── scan_h2lmm.R                         # Kept (not replaced)
└── archive_replaced_by_packrat/         # Old helpers (archived)
    ├── plots.R
    └── tables.R
```

### Key Changes in `21_makeLociPackets.R`

1. **Loads PackRat modules** at top of script (lines 25-30)
2. **Removed inline functions** that are now in PackRat:
   - `get_genes_in_region()` → `join_gene_tables.R`
   - `extract_human_cardiac_traits()` → `phenotype_extract.R`
   - `extract_mouse_cardiac_traits()` → `phenotype_extract.R`
   - `format_human_only()`, `format_mouse_only()`, `format_both()` → `phenotype_extract.R`
   - `generate_gene_info_excel()` → `excel_generation.R`
   - `create_locus_readme()` → `summary_reports.R`

3. **Uses PackRat functions**:
   - `cardiac_keywords()` - predefined phenotype keywords
   - `build_gene_table()` - joins relational tables
   - `add_locus_membership()` - adds locus presence columns
   - `reorder_gene_columns()` - consistent column ordering
   - `create_gene_workbook()` - generates Excel files
   - `create_variant_table()` - generates SNP tables
   - `extract_human_traits()`, `extract_mouse_phenotypes()` - phenotype filtering
   - `create_locus_id()` - file-safe locus identifiers

4. **Plotting**: Currently still uses inline code (marked with TODO for future refactoring)

### What's Preserved

- ✅ All CC-specific data loading logic
- ✅ QTL scan data structures
- ✅ Locus clustering algorithm
- ✅ Expression data integration (VST by Drug/Sex)
- ✅ Founder variant tables
- ✅ Main loop structure
- ✅ Final zip file creation
- ✅ All command-line options

### What's Improved

- ✅ ~480 lines removed from main script
- ✅ Reusable functions in separate modules
- ✅ Better documentation with roxygen comments
- ✅ Easier to test individual components
- ✅ Generalizable to other studies (see `scripts/packrat/README.md`)

## Testing on Slurm Cluster

### Prerequisites

```bash
# Activate conda environment
set +u
source ~/mambaforge/etc/profile.d/conda.sh
conda activate miqtl-env
set -u
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib"
```

### Quick Test (Subset Mode)

Test with a small subset to verify the refactoring works:

```bash
cd /path/to/cc_gwas

# Run on just one locus cluster
sbatch scripts/packets/snakemake/snake.sh
```

**What to check:**
1. ✅ Script loads without errors
2. ✅ PackRat modules load successfully
3. ✅ Excel files generated with correct structure
4. ✅ READMEs created with cardiac gene summaries
5. ✅ Plots generated (if applicable)
6. ✅ Final zip file created

**Expected output:**
```
results/qtl_packets/
├── locus_chr3_131-150Mb/
│   ├── README_summary.txt
│   ├── gene_info_cluster_chr3_131-150Mb.xlsx
│   ├── cardiac_genes_summary.xlsx
│   ├── founder_snp_table_chr3_131-150Mb.csv
│   └── zoomPlots/
│       └── locus_zoom_*.pdf
└── all_loci.zip
```

### Compare Outputs

To verify the refactored version produces identical results:

```bash
# Run original version (if needed for comparison)
Rscript scripts/packets/21_makeLociPackets_ORIGINAL.R \
  --input_sig_regions data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv \
  --output_zip results/qtl_packets/original_all_loci.zip

# Run refactored version
Rscript scripts/packets/21_makeLociPackets.R \
  --input_sig_regions data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv \
  --output_zip results/qtl_packets/refactored_all_loci.zip

# Compare Excel files (manual inspection)
# Check that gene counts, columns, and annotations match
```

### Full Pipeline Test

Once quick test passes, run full pipeline:

```bash
export SNAKEMAKE_MODE="full"
sbatch scripts/packets/snakemake/snake.sh
```

### Troubleshooting

**If modules fail to load:**
- Check that you're on the correct branch: `git branch`
- Verify `scripts/packrat/` directory exists
- Check R can find the modules: `ls scripts/packrat/*.R`

**If functions not found:**
- Ensure PackRat modules are sourced before use
- Check for typos in function names
- Verify function exists in correct module (see `scripts/packrat/README.md`)

**If Excel files look different:**
- Check column ordering (may differ from original)
- Verify conditional formatting applied to "In_*" columns
- Compare sheet names (consolidated vs individual sheets)

**If errors about missing data:**
- Verify all input files exist and paths are correct
- Check that merged_gene_info table has expected columns
- Ensure CC variant tables are in correct location

## Known Issues / TODOs

1. **Plotting**: Still uses inline code from original script. Need to integrate `plot_locus_zoom()` from `scripts/packrat/plot_locus.R` once tested.

2. **Expression data integration**: Currently hard-coded for CC VST data with Drug/Sex grouping. Could be generalized further.

3. **Cardiac keyword filtering**: Uses CC-specific cardiac terms. Other studies may want different phenotype keywords.

## Next Steps

### For CC Cardiac Study

1. **Test on cluster** with subset of loci
2. **Verify outputs** match original pipeline
3. **Run full pipeline** if test passes
4. **Replace inline plotting** with `plot_locus_zoom()` (optional)

### For Paper Publication

1. **Add example vignette** showing non-CC usage
2. **Create synthetic test data** for PackRat
3. **Add unit tests** for core functions
4. **Submit as GitHub release**

### For Generalization

PackRat is now ready for others to use! See `scripts/packrat/README.md` for:
- Example usage with minimal gene lists
- Adding expression/eQTL data
- Custom phenotype filtering
- Adapting to different organisms

## Questions?

Contact:
- Brian Gural (refactoring)
- Christoph Rau (original pipeline)
- GitHub Issues: https://github.com/RauLabUNC/cc_gwas/issues
