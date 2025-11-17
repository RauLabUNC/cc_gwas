## In this folder, there are 2 subfolders
1.  **barebone** (with its own README file t explain the pipeline)
        *Purpose*: Scripts with fundamental functions for: demultiplexing, transforming gene expression, perform differential expression analysis, perform genome scan, eQTL mapping then generating a list of eQTLs overlapping with cardiac trait QTLs. These output will be used to generate the final loci packages. Users can manually input gene expression & genotype
        *State of the folder*: Minimal convey the backbone of gene expression analysis. The GE_eQTL (parts A and B) scripts will work for individual genes, but not yet automated to scan all genes in parallel.
2.  **snakemake** (in progress)
        *Purpose*: Same principles from **barebone** but this pipeline is for streamlining the process to run all the scripts 