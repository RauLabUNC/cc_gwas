# This script processes pre-loaded scan results and thresholds 
# to find significant regions and generate genome plots.

library(tidyverse)
library(zoo)       # For rollapply
library(miqtl)     # For genome.plotter.whole
library(data.table) # For first() / last()
library(optparse)


option_list <- list(
  make_option(c("--input_scans"), type = "character", help = "Paths to input RDS files of scans"),
  make_option(c("--input_thresholds"), type = "character", help = "Paths to input RDS files of thresholds"),
  make_option(c("--output_summary"), type = "character", help = "Path to output CSV file summarizing significant regions"),
  make_option(c("--output_scans"), type = "character", help = "Path to output RDS file of all scans"),
  make_option(c("--output_thresholds"), type = "character", help = "Path to output RDS file of all thresholds")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)
print(opt)

# --- Load Scan and Threshold Data ---
scans <- lapply(opt$options$input_scans, readRDS)
thresholds <- lapply(opt$options$input_thresholds, readRDS)

drugs <- c("Ctrl", "Iso")

# --- List to Store Significant Region Results ---
all_sig_regions_list <- list()

# --- Iterate Through Traits and Drugs ---

# Combine scan and threshold lists for easier iteration
for (scan in scans) {
  # --- Get Scan and Threshold Data ---
  scan <- scan
  threshold <- thresholds[[which(scans == scan)]] |> as.numeric() # This should be the threshold value(s)
  # --- 1. Find Significant Regions  ---
  # Mark significant markers
  is_sig_marker <- scan$LOD > threshold
  # Group contiguous significant markers
  sig_groups_vector <- rep(NA, length(is_sig_marker))
  tf_b <- FALSE
  group <- 0
  for(i in 1:length(is_sig_marker)){
    tf_a <- tf_b
    tf_b <- is_sig_marker[i]
    if (!is.na(tf_b) && tf_a != tf_b && tf_b == TRUE) { # Added !is.na check
      group <- group + 1
      sig_groups_vector[i] <- group
    } else if (!is.na(tf_b) && tf_b == TRUE) {
      sig_groups_vector[i] <- group
    }
  }
  
  # Create data frame for processing
  sig_df <- data.frame(
    loci = names(scan$LOD),
    LOD = scan$LOD,
    block = sig_groups_vector,
    pos = scan$pos$Mb,
    chr = scan$chr,
    stringsAsFactors = FALSE
  )
  
  # Filter out non-significant markers for block processing
  sig_blocks_df <- sig_df %>% filter(!is.na(block))
  
  if (nrow(sig_blocks_df) > 0) {
    # Calculate peak properties within each block
    peak_info <- sig_blocks_df %>%
      group_by(block) %>%
      summarise(
        chr = data.table::first(chr),
        peak_pos = pos[which.max(LOD)],
        Peak_SNP_ID = loci[which.max(LOD)],
        lead_strain = names(which.max(abs(scan$allele.effects[,Peak_SNP_ID]))),
        max_lod = max(LOD, na.rm = TRUE),
        start_sig_pos = min(pos, na.rm = TRUE), # Start of significant block
        end_sig_pos = max(pos, na.rm = TRUE),   # End of significant block
        .groups = 'drop'
      ) %>%
      mutate(threshold_1.5 = max_lod - 1.5)
    
    # Determine LOD drop bounds using rolling average
    bounds_list <- list()
    for(i in 1:nrow(peak_info)){
      current_block <- peak_info$block[i]
      current_chr <- peak_info$chr[i]
      current_threshold <- peak_info$threshold_1.5[i]
      current_start_sig <- peak_info$start_sig_pos[i]
      current_end_sig <- peak_info$end_sig_pos[i]
      
      # Data for the current chromosome
      chr_df <- sig_df %>% filter(chr == current_chr) %>% arrange(pos)
      
      # Find upper bound (scan upstream from start of sig block)
      upper_candidates <- chr_df %>% 
        filter(pos <= current_start_sig, LOD < current_threshold)
      upper_bound_pos <- if(nrow(upper_candidates) > 0) max(upper_candidates$pos, na.rm=TRUE) else min(chr_df$pos, na.rm=TRUE) # Default to chr start if no marker below threshold found
      
      # Find lower bound (scan downstream from end of sig block)
      lower_candidates <- chr_df %>% 
        filter(pos >= current_end_sig, LOD < current_threshold)
      lower_bound_pos <- if(nrow(lower_candidates) > 0) min(lower_candidates$pos, na.rm=TRUE) else max(chr_df$pos, na.rm=TRUE) # Default to chr end if no marker below threshold found
      
      bounds_list[[i]] <- tibble(
        block = current_block,
        lower_pos_lod_drop = lower_bound_pos,
        upper_pos_lod_drop = upper_bound_pos
      )
    }
    
    bounds_df <- bind_rows(bounds_list)
    
    # Combine peak info with bounds
    sig_regions_summary <- peak_info %>%
      inner_join(bounds_df, by = "block") %>%
      mutate(trait = trait,
              drug = drug_status) %>%
      dplyr::select(trait, drug, block, chr, Peak_SNP_ID, lead_strain,
              upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, 
              max_lod) # Reorder/select final columns
    
    # Add to the master list
    all_sig_regions_list[[paste(trait, drug_status, sep="_")]] <- sig_regions_summary
    
    print(paste("      Found", nrow(sig_regions_summary), "significant region(s)."))
    
  } else {
    print("      No significant regions found above threshold.")
  }
} # End trait loop
# --- Save Combined Significant Regions ---
final_sig_regions_df <- bind_rows(all_sig_regions_list)

# Keep the loci with the widest bounds if multiple overlap
merged_regions <- final_sig_regions_df %>% 
  arrange(trait, drug, chr, upper_pos_lod_drop, desc(lower_pos_lod_drop)) %>%   # sort: leftmost first, widest first
  group_by(trait, drug, chr) %>% 
  mutate(max_end_seen = lag(cummax(lower_pos_lod_drop), default = -Inf),
       keep         = lower_pos_lod_drop > max_end_seen) %>%      # nested = FALSE; outermost = TRUE
  ungroup() %>% 
  filter(keep) %>% 
  dplyr::select(-max_end_seen, -keep)

# --- Save Output ---
# Ensure the output directory exists
summary_dir <- dirname(opt$output_summary)
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
}

# Scans and thresholds are output in same directory
results_dir <- dirname(opt$options$output_scans)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

write.csv(merged_regions, opt$options$output_summary, row.names = FALSE)

saveRDS(boxcox, opt$options$output_scans)
saveRDS(thresh.b, opt$options$output_thresholds)
