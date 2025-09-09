## About
This markdown is here to help track the process of reassembling the QTL mapping pipeline.

## Scripts

#### `SCRIPT_NAME.ext`
*Purpose* Briefly describe what the script does, key options, and when to use it.

*Expected inputs*:
- `path/to/dependency_1` - what it is and why it’s needed
- `path/to/dependency_2` - description
- `path/to/input_{paramA}_{paramB}.ext`- note on templated parameters

*Generated outputs*:
- `path/to/output_{paramA}_{paramB}_{feature}.ext` - description of result
#### `snakemake/01_ropScan.R` 
*Purpose*: Executes a genome scan from the miQTL package on a given trait. Has options for different sample aggregation (i.e. taking the average of a group or not), which treatment group to consider (`drug`), and which normalization to use for the input trait measures (could be removed, since we've already established that Box-Cox is the move)

*Expected Inputs:*
- `scripts/extras/scan_h2lmm.R` - modified form of an `miQTL` function to fix a long forgotten error
- `data/raw/genomes/haplotype_cache_cc_083024` - Happy-style genotype information 
- `data/processed/phenotypes/boxCoxTest/{norm}_{agg}_{drug}.csv` - {brackets} are reflecting the snakemake format for the inputs. This is the phenotype information

*Outputs*: 
- `data/processed/ropscan/{norm}_{agg}_{qtl_trait}_{drug}.rds`

*Action items* 
- Add conditional argument that enables both treatment and control groups to be used when performing scans of baseline traits that were recorded before iso-pump implantation
    - need to specify which traits to use
    - Compare results of each approach (pooled versus split baseline measures)


#### `snakemake/02_permTest.R`
*Purpose* Set significance for genome scans generated in `ropScan.R` by performing many, many scans different random permuations of the trait measures

*Expected inputs*:
- `data/processed/ropscan/{norm}_{agg}_{qtl_trait}_{drug}.rds`
- `data/processed/phenotypes/boxCoxTest/{norm}_{agg}_{drug}.csv`
- `data/raw/genomes/haplotype_cache_cc_083024`

*Outputs*:
- `data/processed/scan_thresholds/{norm}_{agg}_{qtl_trait}_{drug}_threshold.rds` - a single numeric threshold value (not sure if -log10(p-value) or LOD score), above which is significant
- `data/processed/scan_thresholds/{norm}_{agg}_{qtl_trait}_{drug}_perm.rds` - Full permuation object. Currently not again used, but takes a long time to generate and would be useful if we wanted to try new thresholds or compare thresholds between traits


#### `03_detectSigLoci.R`
*Purpose* Integrate outputs of ROP scan and perm test for all traits/drug/norm combinations to define and record regions of the genome significantly associated with a given trait.

*Expected inputs*:
- `data/processed/scan_thresholds/{norm}_{agg}_{qtl_trait}_{drug}_threshold.rds`
- `data/processed/ropscan/{norm}_{agg}_{qtl_trait}_{drug}.rds`

*Outputs*:
- `data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv` - locations and associated traits for all loci identified across all scans
- `results/sig_regions/scan_data.rds` - All scan data
- `results/sig_regions/threshold_data.rds` - All threshold values

*Action items* 
- Check if a rolling average is being used to detect when LOD drops below significance threshold.




#### `04_joinLoci.R`
*Purpose* Formats QTL and eQTL results from miQTL analyses into standardized relational tables. Creates consistent output format for both trait and expression QTL mapping results.

*Expected inputs*:
- `data/processed/joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv` - miQTL trait loci from 03_detectSigLoci.R
- `data/processed/joinLoci/eqtl/miQTL/miQTL_output.csv` - miQTL expression QTL results

*Generated outputs*:
- `data/processed/joinLoci/relational_tables/traitLoci.csv` - Standardized trait QTL table with formatted columns
- `data/processed/joinLoci/relational_tables/expLoci.csv` - Standardized expression QTL table with formatted columns

*Action items*
- Add actual significance thresholds (currently NA) based on permutation testing values
- Extract Peak_SNP_pos_bp for miQTL eQTLs (currently NA) from scan results
- Consider adding gene annotation for eQTLs to identify cis vs trans effects


#### `10_joinPositionID.R`
*Purpose* Consolidates all unique positional IDs from trait and expression QTL tables. Creates a unified position reference table with parsed chromosome coordinates and calculated locus lengths.

*Expected inputs*:
- `data/processed/joinLoci/relational_tables/traitLoci.csv` - Trait QTL table with Position_ID column
- `data/processed/joinLoci/relational_tables/expLoci.csv` - Expression QTL table with Position_ID column

*Generated outputs*:
- `data/processed/joinLoci/relational_tables/pos.csv` - Deduplicated position table with chr, start_bp, end_bp, and length_bp

*Action items*
- Consider adding validation for Position_ID format before parsing


#### `11_mergeLociDetectGenes.R`
*Purpose* Identifies genes within QTL loci using biomaRt gene annotations. Performs genomic overlap analysis between loci positions and mouse gene coordinates, filtering out predicted genes (Gm*, *Rik).

*Expected inputs*:
- `data/processed/joinLoci/relational_tables/pos.csv` - Position reference table from 10_joinPositionID.R
- biomaRt connection to Ensembl mmusculus_gene_ensembl database (requires internet)

*Generated outputs*:
- `data/processed/joinLoci/relational_tables/genesInLoci.rds` - RDS file containing gene lists (symbols and Ensembl IDs) for each locus

*Action items*
- Add error handling for biomaRt connection failures
- Consider caching gene annotations locally for reproducibility
- Add option to include/exclude predicted genes based on analysis needs
- Log statistics about gene overlap counts per locus


#### `21_getMouseGenePhenotypes.R`
*Purpose* Queries MouseMine database for phenotype annotations of genes identified in QTL loci. Uses InterMineR to fetch mammalian phenotype ontology annotations for mouse genes.

*Expected inputs*:
- `results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv` - List of genes with mouse symbols
- MouseMine API connection (requires internet)

*Generated outputs*:
- `data/processed/joinLoci/relational_tables/mouseGenePhenotypes.csv` - Mouse gene phenotype annotations

*Action items*
- Add error handling for API connection failures
- Consider batching queries for large gene lists
- Add option to filter phenotypes by relevance (e.g., cardiac-specific)


#### `20_getPathogenicity.R`
*Purpose* Fetches comprehensive gene information including human orthologs, disease associations, genetic constraints, and drug tractability from Open Targets GraphQL API. Processes genes from loci in batches to build relational tables for downstream analysis.

*Expected inputs*:
- `data/processed/joinLoci/relational_tables/genesInLoci.rds` - Genes per locus from 11_mergeLociDetectGenes.R
- BioMart connection for mouse/human ortholog mapping
- Open Targets GraphQL API connection (requires internet)

*Generated outputs*:
- `data/processed/joinLoci/relational_tables/gene_info.csv` - Basic gene information with orthologs
- `data/processed/joinLoci/relational_tables/associations.csv` - Disease associations from Open Targets
- `data/processed/joinLoci/relational_tables/constraints.csv` - Genetic constraint scores
- `data/processed/joinLoci/relational_tables/tractability.csv` - Drug tractability metrics

*Action items*
- Implement robust error recovery for API failures (partial implementation exists)
- Add progress logging for batch processing
- Consider caching API responses to avoid redundant queries
- Fix batch save logic (currently saves every 20 batches incorrectly)



#### `12_splitGenesToSpecies.R`
*Purpose* Splits the comprehensive gene_info table from getPathogenicity.R into separate species-specific tables for mouse genes and orthology relationships.

*Expected inputs*:
- `data/processed/joinLoci/relational_tables/gene_info.csv` - Combined gene information from 20_getPathogenicity.R

*Generated outputs*:
- `data/processed/joinLoci/relational_tables/genes_mouse.csv` - Mouse-specific gene information with coordinates
- `data/processed/joinLoci/relational_tables/orthology.csv` - Mouse-human orthology mapping table

*Action items*
- Add validation for required columns before splitting
- Consider filtering for quality of orthology relationships


#### `13_multTrait_cis-eQTL_nrvmExp.R`
*Purpose* Creates a filtered gene list for downstream analysis by identifying multi-trait candidate genes with high NRVM expression. Modified version that removes PyLMM cis-eQTL dependencies.

*Expected inputs*:
- `data/processed/joinLoci/relational_tables/genes_mouse.csv` - Mouse gene information
- `data/processed/joinLoci/relational_tables/orthology.csv` - Mouse-human orthology
- `data/processed/joinLoci/relational_tables/associations.csv` - Disease associations
- `data/processed/joinLoci/relational_tables/traitLoci.csv` - Trait QTL locations
- `data/processed/joinLoci/relational_tables/mouseGenePhenotypes.csv` - Mouse phenotypes (optional)
- `data/processed/joinLoci/nrvms/bulk_gene.csv` - NRVM expression counts
- `data/processed/joinLoci/nrvms/phenotypes.csv` - NRVM sample metadata

*Generated outputs*:
- `results/joinLoci/geneTables/multTrait_cis-eQTL_nrvmExp.csv` - Filtered gene list with annotations

*Action items*
- Originally required PyLMM cis-eQTL flags - removed in modified version
- Mouse phenotypes made optional to avoid circular dependency


#### `22_makeLociPackets.R`
*Purpose* Final pipeline step that generates comprehensive QTL packets with locus zoom plots, Excel summaries, and README files for each QTL cluster. Modified to remove PyLMM plotting sections.

*Expected inputs*:
- All relational tables from previous steps
- `results/sig_regions/scan_data.rds` - miQTL scan results
- `results/sig_regions/threshold_data.rds` - Significance thresholds
- `data/processed/joinLoci/bulk_exp/` - Bulk expression data (symlinked)
- `data/processed/joinLoci/nrvms/` - NRVM data (symlinked)

*Generated outputs*:
- `results/QTLpackets/README.md` - Summary of all QTL packets
- `results/QTLpackets/locus_*.pdf` - Locus zoom plots for each QTL
- `results/QTLpackets/locus_*.xlsx` - Excel files with gene information

*Action items*
- PyLMM plotting sections commented out
- May need to adjust plot layout without PyLMM track


## August 29th, 2025 Progress

Major pipeline integration completed:
- Removed all PyLMM dependencies throughout pipeline
- Created modified versions of scripts 13 and 22 without PyLMM
- Fixed circular dependency between scripts 13 and 21
- Integrated makeLociPackets.R as final pipeline step
- Symlinked bulk expression and NRVM data from original project
- Fixed SLURM logging to only output in date-specific directories
- Added LD_LIBRARY_PATH workaround for InterMineR

**Blocking issue**: MouseMine appears to be offline/retired - InterMineR cannot connect for phenotype annotations

## August 27th 
Started this approach. Created the `scripts/snakemake` folder, which is for scripts that were run through snakemake previously.
