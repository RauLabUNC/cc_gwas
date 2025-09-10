library(tidyr)
library(dplyr)
library(sva)
library(textshape)
library(edgeR)
library(stringr)
library(tibble)

#STEP 0 - INPUT DEMULTIPLEX DATA------------------------------------------------
#This file has clean strain names (see demultiplex script in CC-eQTL folder)
Clean_Raw = read.csv("Data/ExpressionData/0_Clean_RawCounts_Info_250311.csv", row.names = 1, check.names = F)

Counts= Clean_Raw[,-c(1:6)]
rownames(Counts)=Clean_Raw$SampleID

#STEP 1 - Keep samples (mouse) with library reads >1M---------------------------

HighReadsMice = Counts %>%
  filter(rowSums(.)>1000000)

write.csv(HighReadsMice, "Data/ExpressionData/1_HighReadsSample_RawCounts_250422.csv")


#Make table of sample info
SampleInfo = Clean_Raw[,c(1:6)]
SampleInfo$Plate = str_extract(SampleInfo$CountRep, "^[^_]+")


#Export the list of kept samples
MiceToKeep = rownames(HighReadsMice)
write.csv(MiceToKeep, "Data/ExpressionData/Samples_ReadDepth-1mil.csv", row.names = F)

#Export the updated sample info list
UpdatedSampleInfo = SampleInfo[which(SampleInfo$SampleID %in% MiceToKeep),]
write.csv(UpdatedSampleInfo, "Data/ExpressionData/SampleInfo_ReadDepth-1mil.csv", row.names = F)


#STEP 2 - COMBAT----------------------------------------------------------------

#Raw count matrix after removing small library mice. Still have 55k genes
count_matrix = read.csv("Data/ExpressionData/1_HighReadsSample_RawCounts_250422.csv", check.names = F, row.names = 1)
count_matrix = as.matrix(t(count_matrix))

#Sample info
SampleInfo = read.csv("Data/ExpressionData/SampleInfo_ReadDepth-1mil.csv")


#Correct batch effect based on *plate batches*
#Run ComBat_seq and get adjusted RAW counts. Genes = rows, samples = cols
Adjusted = ComBat_seq(count_matrix, batch = SampleInfo$Plate, group = SampleInfo$Drug)

#Do quick check that the samples match
colnames(Adjusted) == SampleInfo$SampleID

write.csv(Adjusted, "Data/ExpressionData/2_ComBatSeq_250422.csv")


#STEP 3 - FILTER FOR GENES WITH CPM >1 = ~13k genes-----------------------------
#

Adjusted = read.csv("Data/ExpressionData/2_ComBatSeq_250422.csv", check.names = F, row.names = 1)

#CPM, row = genes, col = samples
DGE = DGEList(counts=Adjusted)
CPM = cpm(DGE, log = FALSE, normalized.lib.sizes=TRUE)

#Keep genes with average CPM>1
HighCPM = as.data.frame(CPM) %>%
  filter(rowMeans(.)>1)
NewGeneList = rownames(HighCPM)

write.csv(NewGeneList, "Data/ExpressionData/3_13kGeneNames.csv")

Filter_adj = Adjusted[which(rownames(Adjusted) %in% NewGeneList),]
write.csv(Filter_adj, "Data/ExpressionData/3_13k_RawCounts.csv")


#STEP 4 - TABLE OF UNTRANSFORMED FILTERED DATA----------------------------------

#Sample info
SampleInfo = read.csv("Data/ExpressionData/SampleInfo_ReadDepth-1mil.csv")

Counts_Info = Filter_adj %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(x=., y=SampleInfo, join_by(SampleID)) %>%
  select(SampleID, Strain, Sex, Drug, CountRep, Batch, Plate, everything())
write.csv(Counts_Info, "Data/ExpressionData/4_Untransformed_13k_count_info.csv")


#STEP 5 - TRANSFORM DATA--------------------------------------------------------
#Make sure genes = cols for final output

Filter_adj = read.csv("Data/ExpressionData/3_13k_RawCounts.csv", check.names = F, row.names = 1)
SampleInfo = read.csv("Data/ExpressionData/SampleInfo_ReadDepth-1mil.csv")
NewGeneList = read.csv("Data/ExpressionData/3_13kGeneNames.csv")

#MAKE DESEQ OBJECT FROM RAW COUNTS WITH 13k GENES, RUN CTRL VS ISO SEPARATELY FOR VST BUT STILL KEEP SEX AS A COVAR

#Prep DESeqDataSet to input VST (this is not the same one I used for DESeq, that one I merged ctrl & iso while these are separate)

#CTRL ONLY
#count matrix = genes as rows, samples as cols
count_matrix = Filter_adj[,grep("-C",colnames(Filter_adj))]
SampleInfo_group = SampleInfo %>% filter(Drug=="Ctrl")
#SampleInfo_group$SampleID == colnames(count_matrix) #check if smaple IDs match

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
write.csv(vst_merge, "Data/ExpressionData/5d_sepVST_250429.csv") #sepVST means I transformed ctrl and iso separately, this merged table is just to run miQTL more easily

vst_info = vst_merge %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(x=., y=SampleInfo, join_by(SampleID)) %>%
  dplyr::select(SampleID, Strain, Drug, Sex, CountRep, Batch, Plate, NewGeneList$x)  
write.csv(vst_info, "Data/ExpressionData/5d_sepVST_Info_250429.csv") #This is the file with transformed gene expression and smaple info that I eventually used for miQTL



