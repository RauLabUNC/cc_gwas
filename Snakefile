# Define the input and output files
PHENOS = ["BW.day.0", "BW.day.28", "TH", "LV", "RV", "LA", "RA", "Lung", "Liver", "Adrenal", "THW.by.BW.0", "LVW.by.BW.0", "RVW.by.BW.0", "LAW.by.BW.0", "RAW.by.BW.0", "LuW.by.BW.0", "LiW.by.BW.0", "AdrW.by.BW.0", "IVSd.0", "LVIDd.0", "LVPWd.0", "IVSs.0", "LVIDs.0", "LVPWs.0", "HR.0", "EF.0", "FS.0", "LV.Mass.0", "LV.Mass.Corrected.0", "LV.Vold.0", "LV.Vols.0", "LVIDd.by.LVIDs.0", "LVIDd.by.LVIDs.28", "IVSd.28", "LVIDd.28", "LVPWd.28", "IVSs.28", "LVIDs.28", "LVPWs.28", "HR.28", "EF.28", "FS.28", "LV.Mass.28", "LV.Mass.Corrected.28", "LV.Vold.28", "LV.Vols.28", "delta.EF"]

# Define the final target
rule all:
    input:
        expand("data/processed/scans/{pheno}_scan_results.rds", pheno=PHENOS)

# Rule to run miQTL_ROP_scan.R for each phenotype
rule run_miQTL_ROP_scan:
    output:
        "data/processed/scans/{PHENOS}_scan_results.rds"
    resources:
        mem_gb=4
    shell:
        "Rscript scripts/01_01_miQTL_ROP_scan.R {wildcards.PHENOS}"