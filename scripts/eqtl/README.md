# In this folder, there are 3 subfolders

## 1.  **barebone** (with its own README file t explain the pipeline)
-   *Purpose*: Contain scripts with fundamental functions for: demultiplexing, transforming gene expression, perform differential expression analysis, perform genome scan, eQTL mapping then generating a list of eQTLs overlapping with cardiac trait QTLs. These output will be used to generate the final loci packages for gene prioritization/validation. Users can manually input gene expression & genotype
-   *State of the folder*: Minimal convey the backbone of gene expression analysis. The GE_eQTL (parts A and B) scripts will work for individual genes, but not yet automated to scan all genes in parallel.

## 2.  **snakemake** (in progress)
-   *Purpose*: Same principles from **barebone** but this pipeline is for streamlining the process to run all the scripts in order or in parallel given a gene list using a snakemake workflow 
-   *State of the folder*: Pipeline in progress. There are still some absolute folder paths in there as well as manual steps in between scripts

## 2.  **misc**
-   *Purpose*: Miscellaneous old scripts, will delete once the snakemake pipeline is finished