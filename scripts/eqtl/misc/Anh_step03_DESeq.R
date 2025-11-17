library(DESeq2)

#Load the filtered raw counts from step02 (only ~13k genes, batch corrected)
Filter_adj = read.csv("Data/ExpressionData/3_13k_RawCounts.csv", check.names = F, row.names = 1)
count_matrix = Filter_adj

#Put this in a DESeq count matrix
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = SampleInfo, 
  design = ~Sex+Drug
)
#Save this DESeq object just in case
save(dds, file="Data/ExpressionData/DESeq2_dds_250422.RData")
#Can load the object back in if R crashes
load("Data/Processed/DESeq2_dds_250422.RData")

# Run DESeq 
dds.run<- DESeq(dds)

#View results
resultsNames(dds.run)

log2Fold <- results(dds.run, name = "Drug_Iso_vs_Ctrl")[["log2FoldChange"]] #log2foldchange
p_val <- results(dds.run, name = "Drug_Iso_vs_Ctrl")[["pvalue"]] #padj(p value adjusted)

L2FC_df = data.frame(Gene=rownames(results(dds.run, name = "Drug_Iso_vs_Ctrl")), log2Fold, p_val)
write.csv(L2FC_df, "Results/DESeq2/deseq2_13kGenes_250521.csv") #this is the file I use for my loci packet