# --- Load Libraries, data, and args ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(optparse)
  library(sva)
  library(edgeR)
  library(DESeq2)
})

# Define the command-line arguments
option_list <- list(
  make_option(c("--output_vst_counts"), type = "character", 
  help = "Path to output processed counts file", metavar = "FILE"))

# Parse the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Read the input phenotype file
raw_counts <- read.csv("data/processed/expression/cc_raw_counts.csv", row.names = 1, check.names = FALSE)

sampleInfoBasic <- read.csv("data/processed/phenotypes/sampleInfoBasic.csv", stringsAsFactors = FALSE)

# --- Filter to one replicate per sample based on read depth ---

#Make a column that matches count colnames
sampleInfo = sampleInfoBasic %>%
  mutate(CountRep = paste(PlateNumber, AlitheaBarcode, sep="_"),
         gwas_temp_id = paste0(Strain, Sex, Drug))
# Total reads per replicate (column) without transposing the whole matrix
total_reads <- colSums(raw_counts)              # names = CountRep

rep_depth_df <- tibble(
  CountRep   = names(total_reads),
  TotalReads = as.numeric(total_reads)
)

# Join to sample metadata and choose the top-depth replicate per SampleID
rep_choice <- sampleInfo %>%
  inner_join(rep_depth_df, by = "CountRep") %>%
  group_by(SampleID) %>%
  slice_max(order_by = TotalReads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(TotalReads >= 1e6) # filter out samples with < 1 million reads

# Subset count matrix to the chosen replicates
keep <- colnames(raw_counts) %in% rep_choice$CountRep # sanity check, should be all TRUE
filtered_counts <- as.data.frame(raw_counts)[, keep]

# Reformat to sample rows + metadata
counts_long <- as.data.frame(t(filtered_counts))
counts_long$CountRep <- rownames(counts_long)

counts_info <- rep_choice %>%
  dplyr::select(SampleID, Strain, Sex, Drug, CountRep, PlateNumber, TotalReads) %>%
  inner_join(counts_long, by = "CountRep") %>%
  relocate(SampleID, Strain, Sex, Drug, CountRep, PlateNumber, TotalReads)

# Rename counts matrix with sampleID matching the current CountRep
colnames(filtered_counts) <- counts_info$SampleID[match(colnames(filtered_counts), counts_info$CountRep)]

#### --- ComBat normalization --- ####
cb_counts <- sva::ComBat_seq(filtered_counts, batch = counts_info$PlateNumber, group = counts_info$Drug)

#all(colnames(cb_counts) == counts_info$CountRep) # sanity check, should be all TRUE

# --- Filter for genes with average CPM > 1 ---
dge <- edgeR::DGEList(counts = cb_counts)
cpm <- edgeR::cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)

# select genes with average CPM > 1
high_cpm <- as.data.frame(cpm) %>%
  filter(rowMeans(.) > 1)
new_gene_list <- rownames(high_cpm) # ~13k genes
filter_cb <- cb_counts[which(rownames(cb_counts) %in% new_gene_list), ]

# --- VST Norm ---
## Control first
control_cb <- filter_cb[which(counts_info$Drug == "Ctrl"), ]
sample_info_ctrl <- counts_info %>%
  filter(Drug == "Ctrl") %>%
  arrange(SampleID) %>%
  mutate(Sex = as.factor)

# Ensure columns are in the same order
control_cb <- control_cb[, sample_info_ctrl$CountRep]
all(colnames(control_cb) == sample_info_ctrl$CountRep) # sanity check, should be all TRUE

dds <- DESeqDataSetFromMatrix(
  countData = control_cb,
  colData   = sample_info_ctrl,
  design = ~ Sex
)
#Input count matrix = genes as rows, samples as cols
vst <- varianceStabilizingTransformation(dds, blind = F, fitType = "parametric")

vst_ctrl <- vst@assays@data@listData[[1]] %>%
  t() %>%
  as.data.frame()

## Isoproterenol next
iso_cb <- filter_cb[which(counts_info$Drug == "Iso"), ]
sample_info_iso <- counts_info %>%
  filter(Drug == "Iso") %>%
  arrange(SampleID) %>%
  mutate(Sex = as.factor)

# Ensure columns are in the same order
iso_cb <- iso_cb[, sample_info_iso$CountRep]
all(colnames(iso_cb) == sample_info_iso$CountRep) 

dds <- DESeqDataSetFromMatrix(
  countData = iso_cb,
  colData   = sample_info_iso,
  design = ~ Sex
)

vst <- varianceStabilizingTransformation(dds, blind = F, fitType = "parametric")

vst_iso <- vst@assays@data@listData[[1]] %>%
  t() %>%
  as.data.frame()


#### JUST NEED TO ADD THE LAST MERGE PART ####
#### AND GET THE LIBRARIES INSTALLED/LOADED ####


# --- Save Output ---
# Extract output file path from command-line arguments
output_file <- opt$options$output_boxcox

# Ensure output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  print(paste("Creating output directory:", output_dir))
  dir.create(output_dir, recursive = TRUE)
}

print(paste("Saving processed data to:", output_file))
write.csv(pheno_processed, file = output_file, row.names = F)