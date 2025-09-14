#The whole miQTL, permutation, haplotype block pipeline

#The output here is the significant haplotype table containing blocks that are eQTL hits for 1 gene

#This part is snakemake heavy, please see Anh_step04_miqtl_snakemake folder
  #Scripts order 01_01, 02_01c, 03_01b, 03_02. Then use this output for 04_02d, 04_03 (step05)
  #Scripts 03_02 technically isn't part of snakemake workflow, I just run it last to make a giant table with all genes' SNPs in there for the next step

#Main changes from your trait QTL code
  #Linear model when running genome scan  
~ 1 + Sex (no Drug here bc I ran them separately)

  #Permutation, I did it for all 614 trait loci genes
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan,
                                                      subsample.chr = gene_chromosome, #this is the chromosome that gene is on
                                                      method = "permutation", num.samples = 100,
                                                      use.BLUP = T, model.type = "null")

permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1,
                                      chr = gene_chromosome) #this is the chromosome that gene is on