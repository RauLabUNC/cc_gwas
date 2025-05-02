# This script processes pre-loaded scan results and thresholds 
# to find significant regions and generate genome plots.

library(tidyverse)
library(zoo)       # For rollapply
library(miqtl)     # For genome.plotter.whole
library(data.table) # For first() / last()

traits <- c("BW.day.28", "HR.28", "Lung", "THW.by.BW.0","IVSd.28", "LV.Mass.Corrected.28" , "LV.Mass.Corrected.0")
drugs <- c("Ctrl", "Iso")

boxcox <- c()
zscore <- c()
for(trait in traits){
  for(drug in drugs){
    boxcox[[trait]][[drug]] <- readRDS(paste0("data/processed/ropscan/boxcox_individual_", 
                                              trait, "_", drug, ".rds"))
    zscore[[trait]][[drug]] <- readRDS(paste0("data/processed/ropscan/zscore_individual_", 
                                              trait, "_", drug, ".rds"))
  }
}

thresh.b <- c()
thresh.z <- c()
for(trait in traits){
  for(drug in drugs){
    perms <- readRDS(paste0("data/processed/scan_thresholds/boxcox_individual_", 
                   trait, "_", drug, "_perm.rds"))
    
    thresh.b[[trait]][[drug]] <- get.gev.thresholds(
                                    threshold.scans = perms, 
                                    percentile = 0.90, use.lod = T)
    
    perms <- readRDS(paste0("data/processed/scan_thresholds/zscore_individual_", 
                            trait, "_", drug, "_perm.rds"))
    
    thresh.z[[trait]][[drug]] <- get.gev.thresholds(
                                    threshold.scans = perms, 
                                    percentile = 0.90, use.lod = T)
  }
}



# --- Output Directories ---
sig_region_output_dir <- "results/sig_regions"
genome_plot_output_dir <- "results/genome_plots"
dir.create(sig_region_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(genome_plot_output_dir, recursive = TRUE, showWarnings = FALSE)

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


# --- Iterate Through Normalizations, Traits, and Drugs ---

# Combine scan and threshold lists for easier iteration
scan_data <- list(boxcox = boxcox, zscore = zscore)
threshold_data <- list(boxcox = thresh.b, zscore = thresh.z)
norm_types <- names(scan_data)

for (norm_method in norm_types) {
  print(paste("--- Processing Normalization:", norm_method, "---"))
  
  current_scans <- scan_data[[norm_method]]
  current_thresholds <- threshold_data[[norm_method]]
  
  for (trait in traits) {
    print(paste("  Processing Trait:", trait))
    
    for (drug_status in drugs) {
      print(paste("    Processing Drug Status:", drug_status))
      
      # --- Get Scan and Threshold Data ---
      scan <- current_scans[[trait]][[drug_status]]
      threshold_values <- current_thresholds[[trait]][[drug_status]] # This should be the threshold value(s)
      
      # Check if scan or threshold data is missing/invalid
      if (is.null(scan) || is.null(threshold_values) || length(scan$LOD) == 0 || !is.numeric(threshold_values) || length(threshold_values) == 0) {
        warning(paste("Skipping:", norm_method, trait, drug_status, "- Missing or invalid scan/threshold data."))
        next 
      }
      
      threshold <- threshold_values[1] # Take the first value as the threshold
      if (!is.numeric(threshold) || is.na(threshold)) {
        warning(paste("Skipping:", norm_method, trait, drug_status, "- Invalid threshold value:", threshold))
        next
      }
      
      print(paste("      Using LOD Threshold:", round(threshold, 3)))
      
      # --- 1. Find Significant Regions (Adapting user's logic) ---
      
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
                 drug = drug_status,
                 norm = norm_method) %>%
          dplyr::select(norm, trait, drug, block, chr, 
                 upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, 
                 max_lod) # Reorder/select final columns
        
        # Add to the master list
        all_sig_regions_list[[paste(norm_method, trait, drug_status, sep="_")]] <- sig_regions_summary
        
        print(paste("      Found", nrow(sig_regions_summary), "significant region(s)."))
        
      } else {
        print("      No significant regions found above threshold.")
      }
      
      # --- 2. Generate Genome Plot ---
      plot_filename <- file.path(genome_plot_output_dir, paste0(norm_method, "_", trait, "_", drug_status, "_scan.png"))
      print(paste("      Generating plot:", plot_filename))
      
      png(file = plot_filename, width = 8, height = 4, units = "in", res = 300)
      tryCatch({
        # Use the threshold value directly
        miqtl::genome.plotter.whole(scan.list = list(ScanResult = scan), 
                                    hard.thresholds = threshold,
                                    use.lod = T) # Pass the numeric threshold
        title(main = paste(norm_method, trait, drug_status)) # Add title
      }, error = function(e) {
        warning(paste("Could not generate plot for", norm_method, trait, drug_status, ":", e$message))
        # Create a blank plot with error message if plotting fails
        plot(1, type="n", axes=FALSE, xlab="", ylab="")
        text(1, 1, paste("Plotting Error for:", trait, drug_status, "\n", e$message), cex=0.8)
      })
      dev.off()
      
    } # End drug loop
  } # End trait loop
} # End norm loop

# --- Save Combined Significant Regions ---
if (length(all_sig_regions_list) > 0) {
  final_sig_regions_df <- bind_rows(all_sig_regions_list)
  output_csv_path <- file.path(sig_region_output_dir, "all_significant_regions_summary.csv")
  print(paste("Saving combined significant regions summary to:", output_csv_path))
  write.csv(final_sig_regions_df, output_csv_path, row.names = FALSE)
  saveRDS(scan_data, file = "results/sig_regions/scan_data.rds")
  saveRDS(threshold_data, file = "results/sig_regions/threshold_data.rds")
} else {
  print("No significant regions found across all analyses to save.")
}