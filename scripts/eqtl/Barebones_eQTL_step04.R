remotes::install_github("gkeele/miqtl")
library(miqtl) 
library(tidyverse)

source("/proj/raulab/projects/cc_gwas/scripts/extras/scan_h2lmm.R") #Brian do we need to do sth about your modified function?

One_Gene = "A gene we got from the trait QTL loci"

#####GENOME SCAN----------------------------------------------------------------
genomecache <- "GENOME_CACHE" #Genome cache file


#Load in phenotype data = gene expression, this includes sample information in the first few columns. Gene names are columns, samples are rows
phenotypes <- EXPRESSION_DATA


#Run genome scan on control group or ISO group
miqtl.rop.scan.scaled <- scan.h2lmm.test(
  genomecache = genomecache,
  data = phenotypes,
  pheno.id="SampleID",
  geno.id="Strain",
  formula = One_Gene ~ 1 + Sex, #One_Gene = expression of a particular gene
  use.multi.impute = F,
  return.allele.effects = T,
  use.fix.par = T,
  chr = "all")

# Save the output as an RDS file
saveRDS(miqtl.rop.scan.scaled, SCAN_FILE)


#####PERMUTATION----------------------------------------------------------------
# Load the genome scan output file
scan_file <- SCAN_FILE 
scan <- readRDS(scan_file)

# Load phenotypes (same as above)
phenotypes <- EXPRESSION_DATA
# Load genome cache (same as above)
genomecache <- GENOME_CACHE


#Get the gene chromosome number for this particular gene
gene_pos = read.csv("A_Reference_Table.csv") %>%
  filter(mouse_gene_symbol == One_Gene)
gene_chromosome <- gene_pos$chr


# Permute phenotypes and run scans on each
permuted_phenotype <- generate.sample.outcomes.matrix(scan.object = scan,
                                                      subsample.chr = gene_chromosome,
                                                      method = "permutation", num.samples = 100,
                                                      use.BLUP = T, model.type = "null")


permuted_scans <- run.threshold.scans(sim.threshold.object = permuted_phenotype, 
                                      keep.full.scans=TRUE,
                                      genomecache  = genomecache, 
                                      data = phenotypes,
                                      use.multi.impute = FALSE, 
                                      scan.seed = 1,
                                      chr = gene_chromosome)

#Obtain the top 85% LOD score. Specifiy use.lod = T if we want LOD score, use.lod = F will produce a p-value number
permute_threshold <- get.gev.thresholds(threshold.scans = permuted_scans, percentile = 0.85, use.lod = T)


# Save the threshold and scan as RDS files
saveRDS(permute_threshold, "output_threshold.RData")
saveRDS(permuted_scans, "output_scan.RData")



#####OBTAIN SNPS, FOUNDER EFFECTS FROM SIGNIFICANT HAPLOTYPE BLOCK--------------
# Load necessary libraries
library(zoo)
library(miqtl)
library(dplyr)


# Load scans/thresholds
treat <- One_Gene
trait <- "Treatment_group"

# Read in scan
scan.file <- GENOME_SCAN
scan <- readRDS(scan.file)

# Read and process threshold data
scan.threshold <- readRDS("output_scan.RData")
threshold <- get.gev.thresholds(scan.threshold, use.lod = TRUE, percentile = 0.85, type = "min")

scan$sig <-  scan$LOD > threshold


#### Define loci bounds (LOD 1.5 difference) ####

# Group each continuous block of significant loci
# initialize a list and initial values
sig.groups <- rep(NA, length(scan$sig))
tf.b <- FALSE
group <- 0
for(i in 1:length(scan$sig)){
  tf.a <- tf.b
  tf.b <- scan$sig[i]
  if(tf.a != tf.b & tf.b == T){
    group <- group + 1
    sig.groups[i] <- group
  }else(
    if(tf.b == T){
      sig.groups[i] <- group
    }
  )
}

# Make df of LOD, group, and loci names
sig.df <- data.frame(loci = names(scan$LOD),
                     LOD = scan$LOD,
                     block = sig.groups,
                     pos = scan$pos$Mb,
                     chr = scan$chr) %>%
  na.omit() %>%
  mutate(chr_block = paste(chr, block, sep="_")) 
sig.df$block <- as.numeric(factor(sig.df$chr_block))

loci.order <- sig.df$loci

# For each group, find the max/min
# Get cutoffs from whichever is lower: max - 1.5 or max - min

sig.df <- sig.df |> 
  group_by(block) |> 
  mutate(max = max(LOD),
         upstream.min = data.table::first(LOD),
         downstream.min = data.table::last(LOD)) |> 
  mutate(scan.up = case_when(
    max - upstream.min < 1.5 ~ T,
    .default = F),
    scan.down = case_when(
      max - downstream.min < 1.5 ~ T,
      .default = F),
    threshold = max - 1.5) |> 
  right_join(sig.df) |> 
  as.data.frame()

row.names(sig.df) <- sig.df$loci

# reorder df 
sig.df <- sig.df[loci.order, ]


#### Use a rolling average of LOD scores


# Define the rolling average function
rolling_avg <- function(x, n = 5) {
  rollapply(x, width = 2*n + 1, FUN = mean, fill = NA, align = "center")
}

# Apply the rolling average function to the LOD column
sig.df$rolling_avg_LOD <- rolling_avg(sig.df$LOD, n = 5)

# Reference scan.up/down to selectively run through each LOD value until..
# .. a value is found that is <= the threshold value
groups <- sig.df |> 
  filter(!is.na(block)) |> 
  group_by(block) |> 
  slice_head() |> 
  dplyr::select(block, chr, threshold)

groups$upper <- groups$lower <- groups$upper.pos <- groups$lower.pos <- NA


for(i in groups$block){
  scan.directions <- sig.df[which(sig.df$block == i), c("scan.up", "scan.down")][1,] |> unlist()
  if(scan.directions[1] == T){
    # Check if the block is at the start of the chr
    upstream.pos <- sig.df |> filter(block == i) |> arrange(pos) |> slice_head() |> pull(pos)
    loci.chr <- sig.df |> filter(block == i) |> slice_head() |> pull(chr)
    chr.start <- sig.df |> filter(chr == loci.chr) |> arrange(pos) |> slice_head() |> pull(pos) 
    if(upstream.pos == chr.start){
      # assign the significant upper bound as the upper bound for the confidence interval
      upper.bound <-  sig.df |> filter(block == i) |>  arrange(pos) |> slice_head()
      groups[i, "upper"] <- upper.bound$loci
      groups[i, "upper.pos"] <- upper.bound$pos
    }else{
      print("we'll scan up")
      # Find the closest upstream position that meets the threshold
      upper.bound <- sig.df |> filter(
        chr == groups$chr[[i]] & pos < upstream.pos & 
          rolling_avg_LOD < groups$threshold[[i]]) |> 
        arrange(pos) |> 
        slice_tail() 
      if(nrow(upper.bound) == 0) {
        upper.bound <- sig.df |> 
          filter(chr == groups$chr[[i]] & pos < upstream.pos) |> 
          arrange(pos) |> 
          slice_head()
      }
      groups[i, "upper"] <- upper.bound$loci
      groups[i, "upper.pos"] <- upper.bound$pos
    }
  }else{
    # assign the significant upper bound as the upper bound for the confidence interval
    upper.bound <-  sig.df |> filter(block == i) |>  arrange(pos) |> slice_head()
    groups[i, "upper"] <- upper.bound$loci
    groups[i, "upper.pos"] <- upper.bound$pos}
  
  ## Now check the downstream positions
  if(scan.directions[2] == T){
    # Check if the block is at the end of the chr
    downstream.pos <- sig.df |> filter(block == i) |> arrange(pos) |> slice_tail() |> pull(pos)
    loci.chr <- sig.df |> filter(block == i) |> slice_tail() |> pull(chr) 
    chr.end <- sig.df |> filter(chr == loci.chr) |> arrange(pos) |> slice_tail() |> pull(pos) 
    if(downstream.pos == chr.end){
      # assign the significant lower bound as the lower bound for the confidence interval
      lower.bound <-  sig.df |> filter(block == i) |> arrange(pos) |> slice_tail()
      groups[i, "lower"] <- lower.bound$loci
      groups[i, "lower.pos"] <- lower.bound$pos
    }else{
      # Find the closest downstream position that meets the threshold
      print("we'll scan down")
      lower.bound <- sig.df |> filter(
        chr == groups$chr[[i]] & pos > downstream.pos & 
          rolling_avg_LOD < groups$threshold[[i]]) |> 
        arrange(pos) |> 
        slice_head() 
      # Add ANOTHER rule to make the lower bound = the last loci on the chr if lower.bound returns no rows
      if(nrow(lower.bound) == 0) {
        lower.bound <- sig.df |> 
          filter(chr == groups$chr[[i]] & pos > downstream.pos) |> 
          arrange(pos) |> 
          slice_tail()
      }
      groups[i, "lower"] <- lower.bound$loci
      groups[i, "lower.pos"] <- lower.bound$pos
    }
  }else{
    lower.bound <-  sig.df |> filter(block == i) |> arrange(pos) |> slice_tail()
    groups[i, "lower"] <- lower.bound$loci
    groups[i, "lower.pos"] <- lower.bound$pos}
}

#####
# prep for saving, remove unneeded info
sig.df$is.sig <- sig.df$block > threshold
sig.df <- sig.df |> dplyr::select(-max:-threshold) |> left_join(groups[,c("block", "chr", "lower", "upper", "lower.pos", "upper.pos")])


# For each SNP, get: block index, chr, start and end of that block, LOD of the SNP, trait (gene), treatment

final.df <- sig.df |> 
  filter(!is.na(block)) |>
  group_by(chr_block) |> 
  arrange(desc(LOD)) |> 
  mutate(start = upper.pos,
         end = lower.pos,
         start.snp = upper,
         end.snp = lower,
         trait = trait,
         treatment = treat,
         SNP = loci) |> 
  select(SNP, pos, chr_block, start, end, chr, start.snp, end.snp, LOD, trait, treatment) 


##Figuring out how to get allele effects on the lead SNP
allele.effects <- scan[["allele.effects"]] |> t() |> as.data.frame()
allele.effects <-  mutate_all(allele.effects, abs)

allele.maxes <- data.frame(snp = row.names(allele.effects) )
allele.maxes$lead.strain <- colnames(allele.effects)[max.col(allele.effects)]


#Get the gene chromosome number
gene_pos = read.csv("A_Reference_Table.csv") %>%
  filter(mouse_gene_symbol == args[1])
gene_chromosome <- as.character(gene_pos$chr)


final.final.df = final.df %>%
  left_join(x=., y=allele.maxes, join_by(SNP==snp)) %>%
  filter(chr == gene_chromosome) %>%
  dplyr::select(-chr_block, -start.snp, -end.snp)

final.final.df$pos = final.final.df$pos*1000000

if (nrow(final.final.df)==0){
  final.sig.df = final.final.df
  final.sig.df[1,"trait"] = args[1]
} else{final.sig.df = final.final.df}


write.csv(final.sig.df, "ONE_GENE_ALL_SNPS_FILE.csv", row.names = F) #There will be many files, 1 file per gene. They will all follow the same name pattern, gene's name is in file name


#####CONSOLIDATE ALL GENES, ALL SNPS TO ONE TABLE-------------------------------
library(dplyr); library(tidyverse)

#Do this separately for control and ISO groups
setwd(file.path("Where you've saved all the the ONE_GENE_ALL_SNPS_FILE.csv files"))
MyFiles <- list.files(pattern="*.csv", all.files = T)
df = do.call(rbind, lapply(MyFiles, read.csv))
write.csv(df, file="control_snps.csv") #This can also be iso_snps.csv

