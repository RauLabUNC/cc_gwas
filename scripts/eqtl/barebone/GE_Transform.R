library(tidyr)
library(dplyr)
library(sva)
library(textshape)
library(edgeR)
library(stringr)
library(tibble)

#Note: this long script contains all the steps to process the raw count data from demultiplexing to final VST transformed data with sample info for subsequent eQTL mapping. Each step will produce a csv file that can be used as input for the next step. R Studio might occassionally crash between steps, but we can load in the file from the previous step and continue. 
 
#STEP 0 - INPUT DEMULTIPLEXED DATA----------------------------------------------
#This file has sample ID as rows, gene expression as columns, plus some sample information columns (following demultiplex script)
Clean_Raw = read.csv("Demultiplexed_Raw_Count.csv", row.names = 1, check.names = F)

Counts= Clean_Raw[,-c(1:6)] #Remove sample info columns, keep only gene expression columns
rownames(Counts)=Clean_Raw$SampleID

#STEP 1 - Keep samples (mouse) with library reads >1M---------------------------
HighReadsMice = Counts %>%
  filter(rowSums(.)>1000000)

write.csv(HighReadsMice, "1_HighReadsSample_RawCounts.csv")


#Make table of sample info, extract this information form the first few columns in the demultiplexed file
SampleInfo = Clean_Raw[,c(1:6)]
SampleInfo$Plate = str_extract(SampleInfo$CountRep, "^[^_]+")


#Export the list of samples we kept
MiceToKeep = rownames(HighReadsMice)
write.csv(MiceToKeep, "Samples_ReadDepth-1mil.csv", row.names = F)

#Export the updated sample info list
UpdatedSampleInfo = SampleInfo[which(SampleInfo$SampleID %in% MiceToKeep),]
write.csv(UpdatedSampleInfo, "SampleInfo_ReadDepth-1mil.csv", row.names = F)


#STEP 2 - COMBATSEQ-------------------------------------------------------------

#Raw count matrix after removing small library mice. We still have all ~55000 genes
count_matrix = read.csv("1_HighReadsSample_RawCount.csv", check.names = F, row.names = 1)
count_matrix = as.matrix(t(count_matrix)) #make sure the matrix is genes = rows, samples = cols so that ComBat_seq works

#Sample info
SampleInfo = read.csv("SampleInfo_ReadDepth-1mil.csv")

#Correct batch effect based on *plate batches*
#Run ComBat_seq and get adjusted RAW counts. Genes = rows, samples = cols
Adjusted = ComBat_seq(count_matrix, batch = SampleInfo$Plate, group = SampleInfo$Drug)

#Do a quick check that the samples match
colnames(Adjusted) == SampleInfo$SampleID

#Save the adjusted, batch corrected raw counts
write.csv(Adjusted, "2_ComBatSeq.csv")


#STEP 3 - FILTER FOR GENES WITH CPM >1 = ~13k genes-----------------------------
#Input ComBat_seq adjusted raw counts if needed
Adjusted = read.csv("2_ComBatSeq.csv", check.names = F, row.names = 1)

#Convert expression to CPM, input the matrix as row = genes, col = samples
DGE = DGEList(counts=Adjusted)
CPM = cpm(DGE, log = FALSE, normalized.lib.sizes=TRUE)

#Keep genes with average CPM>1
HighCPM = as.data.frame(CPM) %>%
  filter(rowMeans(.)>1)
NewGeneList = rownames(HighCPM)

#Save the list of genes with average CPM>1 if needed
write.csv(NewGeneList, "3_13kGeneNames.csv")

#Filter the ComBat_seq adjusted raw counts to keep only these ~13k genes
Filter_adj = Adjusted[which(rownames(Adjusted) %in% NewGeneList),]
write.csv(Filter_adj, "3_13k_RawCounts.csv")


#STEP 4 - TABLE OF UNTRANSFORMED FILTERED DATA ALONG WITH SAMPLE INFO-----------

#Sample info
SampleInfo = read.csv("SampleInfo_ReadDepth-1mil.csv")

#Use the batch corrected data
Filter_adj = read.csv("3_13k_RawCounts.csv", check.names = F, row.names = 1)

Counts_Info = Filter_adj %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(x=., y=SampleInfo, join_by(SampleID)) %>%
  select(SampleID, Strain, Sex, Drug, CountRep, Batch, Plate, everything())

write.csv(Counts_Info, "4_Untransformed_13k_count_info.csv") #This is the untransformed, batch corrected raw counts with sample info that we can use for DESeq2 differential expression analysis


#STEP 5 - TRANSFORM DATA--------------------------------------------------------
#Make sure genes = cols for final output

Filter_adj = read.csv("3_13k_RawCounts.csv", check.names = F, row.names = 1)
SampleInfo = read.csv("SampleInfo_ReadDepth-1mil.csv")
NewGeneList = read.csv("3_13kGeneNames.csv")

#Prepare a DESeqDataSet to input into the VST transformation function (this is not the same one DESeq2 DEG analysis, that one has all samples, while these are separate between Control vs ISO treated)

#CTRL ONLY
  #count matrix = genes as rows, samples as cols
count_matrix = Filter_adj[,grep("-C",colnames(Filter_adj))] #Our sample IDs have "-C" or "-I" in their names to indicate control or ISO treatment
SampleInfo_group = SampleInfo %>% filter(Drug=="Ctrl")

SampleInfo_group$SampleID == colnames(count_matrix) #this is just to check if sample IDs match

library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = SampleInfo_group, 
  design = ~Sex
)

#Input count matrix = genes as rows, samples as cols
vst = varianceStabilizingTransformation(dds, blind = F, fitType = "parametric")

vst_ctrl = vst@assays@data@listData[[1]] %>%
  t() %>%
  as.data.frame()


#ISO ONLY
#count matrix = genes as rows, samples as cols
count_matrix = Filter_adj[,grep("-I",colnames(Filter_adj))]
SampleInfo_group = SampleInfo %>% filter(Drug=="Iso")
#SampleInfo_group$SampleID == colnames(count_matrix) #check if sample IDs match

library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = SampleInfo_group, 
  design = ~Sex
)

#Input count matrix = genes as rows, samples as cols
vst = varianceStabilizingTransformation(dds, blind = F, fitType = "parametric")

vst_iso = vst@assays@data@listData[[1]] %>%
  t() %>%
  as.data.frame()


#MERGE 2 GROUPS
vst_merge = rbind(vst_ctrl, vst_iso)
write.csv(vst_merge, "5_VST.csv") #this is just the transformed gene expression table (that can be saved and used for other purposes), we still need to add sample info

vst_info = vst_merge %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(x=., y=SampleInfo, join_by(SampleID)) %>%
  dplyr::select(SampleID, Strain, Drug, Sex, CountRep, Batch, Plate, NewGeneList$x)  
write.csv(vst_info, "5_VST_Info.csv") #This is the file with transformed gene expression and sample info that is eventually used for eQTL mapping



