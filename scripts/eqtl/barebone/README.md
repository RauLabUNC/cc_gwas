# README: Gene Expression Data Analysis Pipeline

This README provides an overview of the pipeline for transforming gene expression data and performing downstream analyses, including eQTL mapping and differential expression analysis

## Purpose of the Pipeline

1.  **Demultiplexing (Step 01)**: Demultiplex raw reads for each sample from pooled libaries based on sample-specific barcode
2.  **Gene Expression Transformation (Step 02)**: Preprocess and transform gene expression data for downstream analyses.
3.  **Differential Expression Analysis (Step 03)**: Use DESeq2 to identify differentially expressed genes.
4.  **eQTL Mapping (Steps 04, 05)**: Identify genetic variants associated with gene expression levels and whether those variants overlap with the trait-significant locus the gene resides within

## Steps in the Pipeline

### Step 01: Raw Count Demultiplexing

-   **Script**: `GE_Demultiplex.R`
-   **Purpose**: Obtain sample-specific raw counts from pooled libraries

### Step 02: Transform Gene Expression Data

-   **Script**: `GE_Transform.R`
-   **Purpose**: Preprocess and transform raw gene expression data to prepare it for both eQTL mapping and differential expression analysis.

### Step 03: Differential Expression Analysis with DESeq2

-   **Script**: `GE_DESEq2.R`
-   **Purpose**: Perform differential expression analysis using DESeq2 to identify genes with significant expression changes in response to ISO.

### Steps 04, 05: eQTL Mapping (following step 02, independent of step 03)

-   **Scripts**:
    -   `GE_eQTL_partA.R`
    -   `GE_eQTL_partB.R`
-   **Purpose**: For each gene within the significant loci hits from trait QTL analyses, perform eQTL mapping to find regions on the genome that correlate with that gene's expression. Then determine if any of those regions overlaps with the trait loci.
    -   Step 04: Perform genome scan, permutation testing to identify significant threshold (based on LOD score), identify significant haplotype blocks and obtain a list of SNPs within them as well as founder line contribution to each block
    -   Step 05: Using the SNPs positions from step 04, as the question "Does my gene has an eQTL within its respective trait locus?" Iterate this for all the genes from trait-significant loci
