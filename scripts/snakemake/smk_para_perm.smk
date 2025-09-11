# --- Configuration ---
import os
configfile: "scripts/snakemake/config.yaml"

# --- Set Up Mode-Specific Variables ---
MODE = os.environ.get("SNAKEMAKE_MODE", "full")
LOG_DIR = os.getenv("LOG_DIR", config["paths"]["log_dir"])

# Get mode-specific parameters from the config
MODE_PARAMS = config["params"][MODE]
QTL_TRAIT = MODE_PARAMS["qtl_trait"]
DRUG = MODE_PARAMS["drug"]
CHUNK_SIZE = MODE_PARAMS["chunk_size"]
NUM_PERMS = MODE_PARAMS["num_perms"]

# Calculate chunks dynamically
PERM_CHUNKS = range(1, (NUM_PERMS // CHUNK_SIZE) + 1)

# Define output prefix based on mode
OUTPUT_PREFIX = "test_" if MODE == "test" else ""
OUTPUT_DIR = config["paths"]["output_dir"]

# --- Target Rule ---
rule all:
    input:
        f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/counts_with_meta.csv"

# --- Rule: Run ROP Scan and Generate Permutations ---
rule run_ropscan_and_gen_perms:
    input:
        script=config["paths"]["scripts"]["ropscan"],
        processed_pheno=f"{config['paths']['phenotypes']}/{{drug}}.csv"
    output:
        scan=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}ropscan/{{qtl_trait}}_{{drug}}.rds",
        perms=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}perm_phenos/{{qtl_trait}}_{{drug}}.rds"
    params:
        qtl_trait="{qtl_trait}",
        drug="{drug}",
        mode=MODE
    resources:
        mem_mb=config["resources"]["run_ropscan_and_gen_perms"][MODE]["mem_mb"],
        time=config["resources"]["run_ropscan_and_gen_perms"][MODE]["time"],
        threads=1
    log:
        stdout=os.path.join(LOG_DIR, "jobs", "run_ropscan_and_gen_perms",
                            "{qtl_trait}_{drug}", "stdout.log"),
        stderr=os.path.join(LOG_DIR, "jobs", "run_ropscan_and_gen_perms",
                            "{qtl_trait}_{drug}", "stderr.log"),
        time=  os.path.join(LOG_DIR, "jobs", "run_ropscan_and_gen_perms",
                            "{qtl_trait}_{drug}", "time.txt")
    benchmark:
        os.path.join(LOG_DIR, "bench",
                    "run_ropscan_and_gen_perms.{qtl_trait}_{drug}.tsv")
    shell:
        """
        mkdir -p "$(dirname {output.scan})" "$(dirname {output.perms})" \
                 "$(dirname {log.stdout})" "$(dirname {log.stderr})" \
                 "$(dirname {log.time})" 

        Rscript {input.script} \
            --input {input.processed_pheno} \
            --output_scan {output.scan} \
            --output_perms {output.perms} \
            --qtl_trait {params.qtl_trait} \
            --drug {params.drug} \
            --mode {params.mode} \
            > {log.stdout} 2> {log.stderr}

        echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-NA}}" >> {log.time} || true
        """

# --- Rule: Run Permutation Scan on a Chunk ---
rule run_permutation_chunk:
    input:
        script=config["paths"]["scripts"]["perm_chunk"],
        perms=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}perm_phenos/{{qtl_trait}}_{{drug}}.rds",
        pheno=f"{config['paths']['phenotypes']}/{{drug}}.csv"
    output:
        scan_chunk=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}perm_chunks/{{qtl_trait}}_{{drug}}/{{chunk}}.rds"
    params:
        qtl_trait="{qtl_trait}",
        drug="{drug}",
        chunk_size=CHUNK_SIZE,
        mode=MODE
    resources:
        mem_mb=config["resources"]["run_permutation_chunk"][MODE]["mem_mb"],
        time=config["resources"]["run_permutation_chunk"][MODE]["time"]
    log:
        stdout=os.path.join(LOG_DIR, "jobs", "run_permutation_chunk",
                            "{qtl_trait}_{drug}", "{chunk}", "stdout.log"),
        stderr=os.path.join(LOG_DIR, "jobs", "run_permutation_chunk",
                            "{qtl_trait}_{drug}", "{chunk}", "stderr.log"),
        time=  os.path.join(LOG_DIR, "jobs", "run_permutation_chunk",
                            "{qtl_trait}_{drug}", "{chunk}", "time.txt")
    benchmark:
        os.path.join(LOG_DIR, "bench",
                    "run_permutation_chunk.{qtl_trait}_{drug}_{chunk}.tsv")
    shell:
        """
        mkdir -p "$(dirname {output.scan_chunk})" \
                 "$(dirname {log.stdout})"

        Rscript {input.script} \
            --input_perms {input.perms} \
            --input_pheno {input.pheno} \
            --output_scan_chunk {output.scan_chunk} \
            --qtl_trait {params.qtl_trait} \
            --drug {params.drug} \
            --chunk_index {wildcards.chunk} \
            --chunk_size {params.chunk_size} \
            --mode {params.mode} > {log.stdout} 2> {log.stderr}

        echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-NA}}" >> {log.time} || true
        """

# --- Rule: Aggregate Permutation Chunks and Get Threshold ---
rule aggregate_permutations:
    input:
        script=config["paths"]["scripts"]["perm_agg"],
        chunks=lambda wildcards: expand(f"{OUTPUT_DIR}/{OUTPUT_PREFIX}perm_chunks/{wildcards.qtl_trait}_{wildcards.drug}/{{chunk}}.rds",
                                        chunk=PERM_CHUNKS)
    output:
        threshold=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}scan_thresholds/{{qtl_trait}}_{{drug}}_threshold.rds",
        max_stats=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}scan_thresholds/{{qtl_trait}}_{{drug}}_max_stats.rds"
    params:
        qtl_trait="{qtl_trait}",
        drug="{drug}"
    resources:
        mem_mb=config["resources"]["aggregate_permutations"][MODE]["mem_mb"],
        time=config["resources"]["aggregate_permutations"][MODE]["time"]
    log:
        stdout=os.path.join(LOG_DIR, "jobs", "aggregate_permutations",
                            "{qtl_trait}_{drug}", "stdout.log"),
        stderr=os.path.join(LOG_DIR, "jobs", "aggregate_permutations",
                            "{qtl_trait}_{drug}", "stderr.log"),
        time=os.path.join(LOG_DIR, "jobs", "aggregate_permutations",
                          "{qtl_trait}_{drug}", "time.txt")
    shell:
        """
        mkdir -p "$(dirname {output.threshold})" "$(dirname {log.stdout})"

        Rscript {input.script} \
            --output_threshold {output.threshold} \
            --output_max_stats {output.max_stats} \
            {input.chunks} > {log.stdout} 2> {log.stderr}
        
        echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-NA}}" >> {log.time} || true
        """

# --- Rule: Detect Significant Loci ---
# This rule aggregates all scan results and thresholds to identify significant QTL regions
rule detect_significant_loci:
    input:
        scans=expand(f"{OUTPUT_DIR}/{OUTPUT_PREFIX}ropscan/{{trait}}_{{drug}}.rds",
                     trait=QTL_TRAIT, drug=DRUG),
        thresholds=expand(f"{OUTPUT_DIR}/{OUTPUT_PREFIX}scan_thresholds/{{trait}}_{{drug}}_threshold.rds",
                          trait=QTL_TRAIT, drug=DRUG),
        script=config["paths"]["scripts"]["detect_loci"]
    output:
        summary=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/trait_qtl/miQTL/all_significant_regions_summary.csv",
        relational_summary=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/traitLoci.csv",
        pos_summary=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/pos.csv",
        all_scans=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}trait_qtl/all_scans.rds",
        all_thresholds=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}trait_qtl/all_thresholds.rds"
    params:
        scans_csv=lambda wc, input: ",".join(input.scans),
        thresholds_csv=lambda wc, input: ",".join(input.thresholds)
    log:
        stdout=os.path.join(LOG_DIR, "jobs", "detect_significant_loci", "stdout.log"),
        stderr=os.path.join(LOG_DIR, "jobs", "detect_significant_loci", "stderr.log"),
        time=os.path.join(LOG_DIR, "jobs", "detect_significant_loci", "time.txt")
    resources:
        time=config["resources"]["detect_significant_loci"][MODE]["time"]
    shell:
        """
        mkdir -p "$(dirname {output.summary})" "$(dirname {log.stdout})"
        Rscript {input.script} \
          --input_scans "{params.scans_csv}" \
          --input_thresholds "{params.thresholds_csv}" \
          --output_summary {output.summary} \
          --output_relational_summary {output.relational_summary} \
          --output_pos_summary {output.pos_summary} \
          --output_scans {output.all_scans} \
          --output_thresholds {output.all_thresholds} > {log.stdout} 2> {log.stderr}
        echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-NA}}" >> {log.time} || true
        """

rule query_genes_in_loci:
    input:
        pos_summary=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/pos.csv",
        script=config["paths"]["scripts"]["query_genes"]
    output:
        genes=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/genesInLoci.rds",
    log:
        stdout=os.path.join(LOG_DIR, "jobs", "query_genes_in_loci", "stdout.log"),
        stderr=os.path.join(LOG_DIR, "jobs", "query_genes_in_loci", "stderr.log"),
        time=os.path.join(LOG_DIR, "jobs", "query_genes_in_loci", "time.txt")
    resources:
        time=config["resources"]["query_genes_in_loci"][MODE]["time"]
    shell:
        """
        mkdir -p "$(dirname {output.genes})" "$(dirname {log.stdout})"
        Rscript {input.script} \
          --input_pos_summary "{input.pos_summary}" \
          --output_genes {output.genes} > {log.stdout} 2> {log.stderr}
        echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-NA}}" >> {log.time} || true
        """

rule make_counts_from_star:
    input:
        genes_in_trait_loci=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/genesInLoci.rds",
        script=config["paths"]["scripts"]["make_counts"]
    output:
        raw_counts=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/raw_counts.csv",
        counts_with_meta=f"{OUTPUT_DIR}/{OUTPUT_PREFIX}joinLoci/relational_tables/counts_with_meta.csv"
    log:
        stdout=os.path.join(LOG_DIR, "jobs", "make_counts_from_star", "stdout.log"),
        stderr=os.path.join(LOG_DIR, "jobs", "make_counts_from_star", "stderr.log"),
        time=os.path.join(LOG_DIR, "jobs", "make_counts_from_star", "time.txt")
    resources:
        time=config["resources"]["query_genes_in_loci"][MODE]["time"]
    shell:
        """
        mkdir -p "$(dirname {output.raw_counts})" "$(dirname {log.stdout})"
        Rscript {input.script} \
          --input_genes_in_trait_loci "{input.genes_in_trait_loci}" \
          --output_raw_counts {output.raw_counts} \
          --output_counts_with_meta {output.counts_with_meta} > {log.stdout} 2> {log.stderr}
        echo "SLURM_JOB_ID=${{SLURM_JOB_ID:-NA}}" >> {log.time} || true
        """