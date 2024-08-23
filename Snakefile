# Define the input and output files
VCF_FILE = "input.vcf"
PHENO_FILE = "phenotypes.txt"
OUTPUT_DIR = "results"


# Define the final target
rule all:
    input:
        directory(OUTPUT_DIR)
        
# Rule to preprocess the VCF file
rule preprocess_vcf:
    input:
        vcf = VCF_FILE
    output:
        vcf = "preprocessed.vcf"
    shell:
        "python preprocess_vcf.py {input.vcf} {output.vcf}"

# Rule to process the phenotypes
rule process_phenotypes:
    input:
        pheno = PHENO_FILE
    output:
        pheno = "processed_phenotypes.txt"
    shell:
        "python process_phenotypes.py {input.pheno} {output.pheno}"

# Rule to run GWAS analysis
rule run_gwas:
    input:
        vcf = "preprocessed.vcf",
        pheno = "processed_phenotypes.txt"
    output:
        gwas_results = directory(OUTPUT_DIR)
    shell:
        "python run_gwas.py {input.vcf} {input.pheno} {output.gwas_results}"
