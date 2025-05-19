# This script is meant to define the confidence intervals for significant loci..
# .. in QTL mapping w/ 1.5 LOD drop from the maximum in a given peak

# Load libs
library("tidyverse")
library("data.table")
library(miqtl)

# Read in the first trailing argument for phenotype
args <- commandArgs(trailingOnly = TRUE)
args <- c("CSA", "iso")



# Scan
scan_file <- file.path("data/processed/scans", args[2], paste0(args[1], "_scan_results.rds"))
scan <- readRDS(scan_file)

# Threshold
threshold_file <- file.path("data/processed/scan_thresholds", args[2], paste0(args[1], "_scan.rds"))
scan_threshold <- readRDS(threshold_file)


threshold <- get.gev.thresholds(
  scan_threshold,
  use.lod = T,
  percentile = 0.85,
  type = c("min")
)
# Label significant loci
scan$sig <-  scan$LOD > threshold

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
                     chr = scan$chr)

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
# Load necessary library
library(zoo)

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
          filter(chr == groups$chr[[i]] & pos < downstream.pos) |> 
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


# prep for saving, remove unneeded info

sig.df$is.sig <- sig.df$block > threshold
sig.df <- sig.df |> dplyr::select(-max:-threshold) |> left_join(groups[,c("block", "chr", "lower", "upper", "lower.pos", "upper.pos")])

# Ensure the directory exists
output_dir <- file.path("data/processed/sig_loci/oneFiveLODRanges", args[2])
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(sig.df, paste0(output_dir,"/", args[[1]], ".csv"), row.names = F)

