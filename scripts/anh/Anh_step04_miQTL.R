#The whole miQTL, permutation, haplotype block pipeline

#The output here is the significant haplotype table containing blocks that are eQTL hits for 1 gene

#This is snakemake heavy, I might need to give you a whole separate folder for it

#Main changes from your trait QTL code
  #Linear model when running genome scan  
~ 1 + Sex (no Drug here bc I ran them separately)

  #Permutation, I did it for all 614 trait loci genes
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan,
                                                      subsample.chr = gene_chromosome, #this is the chromosome that gene is on
                                                      method = "permutation", num.samples = 1000,
                                                      use.BLUP = T, model.type = "null")

permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1,
                                      chr = gene_chromosome) #this is the chromosome that gene is on