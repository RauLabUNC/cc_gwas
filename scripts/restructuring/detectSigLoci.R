# This script processes pre-loaded scan results and thresholds 
# to find significant regions and generate genome plots.

library(tidyverse)
library(zoo)       # For rollapply
library(miqtl)     # For genome.plotter.whole
library(data.table) # For first() / last()

drugs <- c("Ctrl", "Iso")

# List all files in the directory ending with .rds
all_files <- list.files( "data/processed/ropscan/", pattern = "\\.rds$")

# Define regex pattern to capture phenotype and condition (Ctrl or Iso)
# looks for pattern: anything_anything_PHENOTYPE_CONDITION.rds
pattern <- "^.*?_.*?_(.*?)_(Ctrl|Iso)\\.rds$"

# Extract phenotype and condition using the pattern
# str_match returns a matrix: [,1] is full match, [,2] is phenotype, [,3] is condition
matches <- str_match(all_files, pattern)

# Get unique phenotypes where condition is Ctrl (column 3 of matches)
traits <- unique(matches[!is.na(matches[, 3]), 2])

boxcox <- c()
for(trait in traits){
  for(drug in drugs){
      boxcox[[trait]][[drug]] <- readRDS(paste0("data/processed/ropscan/boxcox_individual_", 
                                                trait, "_", drug, ".rds"))
  } 
}

thresh.b <- c()
for(trait in traits){
  for(drug in drugs){
    thresh.b[[trait]][[drug]] <- readRDS(
                                  paste0("data/processed/scan_thresholds/boxcox_individual_", 
                                                           trait, "_", drug, "_threshold.rds"))
    
  }
}



# --- Output Directories ---
sig_region_output_dir <- "results/sig_regions"
dir.create(sig_region_output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Master List to Store Significant Region Results ---
all_sig_regions_list <- list()

# --- Define the rolling average function ---
rolling_avg <- function(x, n = 5) {
  # Ensure there are enough non-NA values for the window
  if (sum(!is.na(x)) < (2*n + 1)) {
    warning("Not enough non-NA values for rolling average window size.")
    return(rep(NA_real_, length(x)))
  }
  zoo::rollapply(x, width = 2*n + 1, FUN = mean, na.rm = TRUE, fill = NA, align = "center", partial = TRUE) # Use partial=TRUE to handle edges
}


# --- Iterate Through Traits and Drugs ---

# Combine scan and threshold lists for easier iteration
  for (trait in traits) {
    print(paste("  Processing Trait:", trait))
    for (drug_status in drugs) {
      # --- Get Scan and Threshold Data ---
      scan <- boxcox[[trait]][[drug_status]]
      threshold <- thresh.b[[trait]][[drug_status]] |> as.numeric()# This should be the threshold value(s)
      
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
        
        # Calculate rolling average LOD for the whole scan first
        sig_df$rolling_avg_LOD <- rolling_avg(sig_df$LOD, n = 5) # Adjust n if needed
        
        
        
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
            filter(pos <= current_start_sig, rolling_avg_LOD < current_threshold)
          upper_bound_pos <- if(nrow(upper_candidates) > 0) max(upper_candidates$pos, na.rm=TRUE) else min(chr_df$pos, na.rm=TRUE) # Default to chr start if no marker below threshold found
          
          # Find lower bound (scan downstream from end of sig block)
          lower_candidates <- chr_df %>% 
            filter(pos >= current_end_sig, rolling_avg_LOD < current_threshold)
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
      
    } # End drug loop
  } # End trait loop
# --- Save Combined Significant Regions ---
final_sig_regions_df <- bind_rows(all_sig_regions_list)
output_csv_path <- file.path("data/processed/joinLoci/trait_qtl/miQTL/", "all_significant_regions_summary.csv")
write.csv(final_sig_regions_df, output_csv_path, row.names = FALSE)
