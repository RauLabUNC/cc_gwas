library(DESeq2)

#Load the filtered raw counts from expression transformation script (with only ~13k genes, batch corrected, VST)
Filter_adj = read.csv("3_13k_RawCounts.csv", check.names = F, row.names = 1)
count_matrix = Filter_adj

#Put this in a DESeq count matrix
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = SampleInfo, 
  design = ~Sex+Drug
)
#Save this DESeq object just in case
save(dds, file="DESeq2_dds.RData")

# Run DESeq2 
dds.run<- DESeq(dds)

#View results
resultsNames(dds.run)

#Extract info on log2 fold change and p value from the DESeq2 results
log2Fold <- results(dds.run, name = "Drug_Iso_vs_Ctrl")[["log2FoldChange"]] #log2 fold change
p_val <- results(dds.run, name = "Drug_Iso_vs_Ctrl")[["pvalue"]] #padj(p value adjusted)

#Make them into a table for easy of reference
L2FC_df = data.frame(Gene=rownames(results(dds.run, name = "Drug_Iso_vs_Ctrl")), log2Fold, p_val)

write.csv(L2FC_df, "deseq2_results.csv") #this is the differential gene expresison file that will be use to assemble gene selection evidence 