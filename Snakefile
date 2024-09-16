# Define the input and output files
PHENOS = ["BW.day.0", "BW.day.28", "TH", "LV", "RV", "LA", "RA", "Lung", "Liver", "Adrenal", "THW.by.BW.0", "LVW.by.BW.0", "RVW.by.BW.0", "LAW.by.BW.0", "RAW.by.BW.0", "LuW.by.BW.0", "LiW.by.BW.0", "AdrW.by.BW.0", "IVSd.0", "LVIDd.0", "LVPWd.0", "IVSs.0", "LVIDs.0", "LVPWs.0", "HR.0", "EF.0", "FS.0", "LV.Mass.0", "LV.Mass.Corrected.0", "LV.Vold.0", "LV.Vols.0", "LVIDd.by.LVIDs.0", "LVIDd.by.LVIDs.28", "IVSd.28", "LVIDd.28", "LVPWd.28", "IVSs.28", "LVIDs.28", "LVPWs.28", "HR.28", "EF.28", "FS.28", "LV.Mass.28", "LV.Mass.Corrected.28", "LV.Vold.28", "LV.Vols.28", "delta.EF", "Wall.Thicknessd.0", "Wall.Thicknessd.28", "Delta.Wall.Thickness.d", "Delta.LVIDs.d", "Delta.FS", "Percent.Fibrosis"]
TREATMENT = ["control", "iso"]
# Define the final  target
rule all:
    input:
        expand("results/genome_scans/{treatment}/{pheno}_scan_results.png", treatment = TREATMENT, pheno=PHENOS),
        expand("results/genome_scans_thresholds/{treatment}/{pheno}_scan_threshold.png", treatment = TREATMENT, pheno=PHENOS),
        expand("data/processed/sig_loci/{treatment}/{pheno}_loci.csv", treatment = TREATMENT, pheno=PHENOS)
# Rule to average phenotypes
rule average_phenotypes:
    output:
        "data/processed/phenotypes/mean_cc_panel_04_16_24.csv"
    shell:
        "Rscript scripts/01_00_average_phenos.R"
# Rule to run miQTL_ROP_scan.R for each phenotype
rule run_miQTL_ROP_scan:
    input:
        "data/processed/phenotypes/mean_cc_panel_04_16_24.csv"
    output:
        "data/processed/scans/{TREATMENT}/{PHENOS}_scan_results.rds"
    resources:
        mem_gb=4
    shell:
        "Rscript scripts/01_01_miQTL_ROP_scan.R {wildcards.PHENOS} {wildcards.TREATMENT}"

rule plot_miQTL_ROP_scan:
    input:
        "data/processed/scans/{TREATMENT}/{PHENOS}_scan_results.rds"
    output:
        "results/genome_scans/{TREATMENT}/{PHENOS}_scan_results.png"
    resources:
        mem_gb=4
    shell:
        "Rscript scripts/01_02_plot_scan.R {wildcards.PHENOS} {wildcards.TREATMENT}"

rule threshold_miQTL_ROP_scan:
    input:
        "data/processed/scans/{TREATMENT}/{PHENOS}_scan_results.rds"
    output:
        "data/processed/scan_thresholds/{TREATMENT}/{PHENOS}_threshold.rds"
    resources:
        mem_gb=4
    shell:
        "Rscript scripts/02_01_permutation.R  {wildcards.PHENOS} {wildcards.TREATMENT}"

rule plot_threshold_scan:
    input:
        "data/processed/scans/{TREATMENT}/{PHENOS}_scan_results.rds",
        "data/processed/scan_thresholds/{TREATMENT}/{PHENOS}_threshold.rds"
    output:
        "results/genome_scans_thresholds/{TREATMENT}/{PHENOS}_scan_threshold.png"
    shell:
        "Rscript scripts/02_02_plot_permutation.R  {wildcards.PHENOS} {wildcards.TREATMENT}"

rule extract_sig_loci:
    input:
        "data/processed/scans/{TREATMENT}/{PHENOS}_scan_results.rds",
        "data/processed/scan_thresholds/{TREATMENT}/{PHENOS}_threshold.rds"
    output:
        "data/processed/sig_loci/{TREATMENT}/{PHENOS}_loci.csv"
    shell:
        "Rscript scripts/03_01_pullSigHaplotypes.R  {wildcards.PHENOS} {wildcards.TREATMENT}"