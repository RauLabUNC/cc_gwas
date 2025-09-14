library(tidyr)
library(dplyr)
library(stringr)
library(edgeR)
library(Matrix)
library(data.table)
library(tidyverse)
library(biomaRt)


################################################################################
#REFORMAT AND MATCH SAMPLE ID
################################################################################

#Set working dir to lab's folder 
setwd("/proj/raulab/projects/alithea_seq")

##Brian's code for de-multiplexing
# Make counts folder
if(!dir.exists("data/processed/counts")){
  dir.create("data/processed/star_aligned/counts_matrices")
}


# Function to read a plate's data
read_plate_data <- function(plate_number) {
  plate_dir <- paste0("data/processed/star_aligned/CC_Plate_", plate_number)
  # Read the sparse matrix
  mat <- readMM(paste0(plate_dir, "/Solo.out/Gene/raw/matrix.mtx"))
  
  # Read features (genes) and barcodes
  features <- fread(paste0(plate_dir, "/Solo.out/Gene/raw/features.tsv"), 
                    header=FALSE, stringsAsFactors=FALSE)
  barcodes <- fread(paste0(plate_dir, "/Solo.out/Gene/raw/barcodes.tsv"), 
                    header=FALSE, stringsAsFactors=FALSE)
  
  # Add plate identifier to barcodes
  barcodes$V1 <- paste0(plate_number, "_", barcodes$V1)
  
  # Convert to dense matrix and add row/column names
  mat <- as.matrix(mat)
  rownames(mat) <- features$V1
  colnames(mat) <- barcodes$V1
  
  return(mat)
}

# Read each plate's data. This is the count matrix, colnames = the plate # _ barcode, rownmes = transcripts
counts <- lapply(1:5, read_plate_data) %>%
  do.call(cbind, .) |> 
  as.data.frame()


#### Adapting the old ens_to_gene ####

# Load the Ensembl dataset for mouse genes
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract the Ensembl gene IDs from the row names
ensembl_ids <- rownames(counts)

# Query biomaRt to get the gene symbols for these Ensembl IDs
genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
               filters = 'ensembl_gene_id', 
               values = ensembl_ids, 
               mart = mart)

# Join the gene symbols back to the counts data frame
# First, make sure the Ensembl IDs are a column in the counts data frame for merging
counts$ensembl_gene_id <- rownames(counts)

# Merge to add gene symbols
counts_with_symbols <- left_join(counts, genes, by = 'ensembl_gene_id')


# First, identify numeric columns once
numeric_cols <- names(counts_with_symbols)[sapply(counts_with_symbols, is.numeric)]

# Convert to data.table
dt <- as.data.table(counts_with_symbols)

# Perform the summarization
summarized_counts <- dt[!is.na(external_gene_name),
                        lapply(.SD, sum, na.rm = TRUE),
                        by = external_gene_name,
                        .SDcols = numeric_cols]

# Replace Ensembl IDs with gene symbols where available, keeping Ensembl ID where not
rownames(summarized_counts) <- summarized_counts$external_gene_name

summarized_counts <- summarized_counts |> dplyr::select(-external_gene_name)

setwd("/proj/raulab/users/Anh/CC-eQTL")
write.csv(summarized_counts, "Data/Processed/ExpressionData/AllRawCounts_250311.csv")



###Anh's part, clean up count matrix more### 

# Import Brian's cleaned up sample info with ID, plate number, barcode, etc
sampleInfoBasic <- read.csv("/proj/raulab/projects/alithea_seq/data/processed/phenos/sampleInfoBasic.csv")
#Make a column that matches count colnames
sampleInfo = sampleInfoBasic %>%
  mutate(CountRep = paste(PlateNumber, AlitheaBarcode, sep="_"))

#Remove technical replicate with lower read depth (total counts) => "clean" counts
counts_2 = as.data.frame(t(summarized_counts)) 
#if R crases, run the line below
counts_2 = read.csv("Data/Processed/ExpressionData/AllRawCounts_250311.csv", check.names = F, row.names = 1)
counts_2 = as.data.frame(t(counts_2))

counts_2$TotalReads = rowSums(counts_2)
counts_2$CountRep = rownames(counts_2)
counts_2 = as.data.frame(counts_2) %>%
  right_join(x=., y=sampleInfo, join_by(CountRep)) %>%
  group_by(SampleID) %>%
  dplyr::slice(which.max(TotalReads))

#The Sample ID aka to-be row names here is by this format: [Strain]-[Drug C or I][Sex F or M][1 or 2 for iso, 1 for ctrl]. Ex: 1-CF1 means strain 1, ctrl, female

#Keep all clean raw counts + condensed sample info
Counts_and_Info = counts_2 %>%
  dplyr::select(-c(TotalReads, PlateNumber, AlitheaBarcode)) %>%
  dplyr::select(SampleID, Strain, Sex, Drug, CountRep, Batch, everything())
write.csv(Counts_and_Info, "Data/Processed/ExpressionData/1_Cleaned_pre-BatchCorrect/Clean_RawCounts_Info_250311.csv") #This file name was changed to "0_Clean_RawCounts_Info_250311" when I moved it to my clean folder