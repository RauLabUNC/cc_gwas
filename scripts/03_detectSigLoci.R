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
  make_option(c("--output_relational_summary"), type = "character", help = "Path to output CSV file for relational table of trait loci"),
  make_option(c("--output_pos_summary"), type = "character", help = "Directory to save position reference table"),
  make_option(c("--output_scans"), type = "character", help = "Path to output RDS file of all scans"),
  make_option(c("--output_thresholds"), type = "character", help = "Path to output RDS file of all thresholds")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, positional_arguments = TRUE)
print(opt)

# --- Load Scans and Thresholds ---
scan_paths <- strsplit(opt$options$input_scans, ",", fixed = TRUE)[[1]]
threshold_paths <- strsplit(opt$options$input_thresholds, ",", fixed = TRUE)[[1]]

if (length(scan_paths) != length(threshold_paths)) {
  stop("Mismatch: number of scans (", length(scan_paths),
       ") != number of thresholds (", length(threshold_paths), ")")
}

scans <- lapply(scan_paths, readRDS)
thresholds <- lapply(threshold_paths, readRDS)

# --- List to Store Significant Region Results ---
all_sig_regions_list <- list()

# --- Iterate Through Traits and Drugs ---
# Combine scan and threshold lists
for (i in seq_along(scans)) {
  scan_path <- scan_paths[i]
  thr_path  <- threshold_paths[i]

  # Derive trait and drug from filenames:
  # ropscan/<trait>_<drug>.rds
  scan_base <- tools::file_path_sans_ext(basename(scan_path))
  parts <- strsplit(scan_base, "_", fixed = TRUE)[[1]]

  trait <- parts[1]
  drug_status <- parts[2]

  print(drug_status)
  print(trait)

  scan <- scans[[i]]
  threshold <- thresholds[[i]]

  is_sig_marker <- scan$LOD > threshold

  sig_groups_vector <- rep(NA, length(is_sig_marker))
  tf_b <- FALSE
  group <- 0
  for (j in seq_along(is_sig_marker)) {
    tf_a <- tf_b
    tf_b <- is_sig_marker[j]
    if (!is.na(tf_b) && tf_a != tf_b && tf_b) {
      group <- group + 1
      sig_groups_vector[j] <- group
    } else if (!is.na(tf_b) && tf_b) {
      sig_groups_vector[j] <- group
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


# Format columns
## Format trait QTL data
trait_loci <- merged_regions |> 
  mutate(Method = "miQTL",
         Analysis_Type = "TraitQTL",
         Dependent_Variable = trait,
         Treatment = drug,
         Chr = chr,
         Locus_Start_bp = upper_pos_lod_drop * 10^6,
         Locus_End_bp = lower_pos_lod_drop * 10^6,
         Peak_SNP_ID = Peak_SNP_ID,
         Peak_SNP_pos_bp = peak_pos * 10^6,
         Lead_Strain = lead_strain,
         Peak_Significance_Value = max_lod,
         Significance_Metric = "LOD",
         Significance_Threshold = NA,
         Locus_ID = paste0(Chr, ":", Locus_Start_bp, "-", Locus_End_bp, "_", 
                           Analysis_Type, "_", Dependent_Variable, "_", Treatment),
         Position_ID = paste0(Chr, ":", Locus_Start_bp, "-", Locus_End_bp)) |>
  dplyr::select(Locus_ID, Position_ID, Method, Analysis_Type, Dependent_Variable, Treatment, Chr,
                Locus_Start_bp, Locus_End_bp, Peak_SNP_ID, Lead_Strain, Peak_SNP_pos_bp, Peak_Significance_Value,
                Significance_Metric, Significance_Threshold)

# Make unique position IDs for relational table
# --- Combine and Parse Position IDs ---
# Union and deduplicate positions
pos <- trait_loci$Position_ID %>% unique() %>%
  as.data.table() %>%
  setnames("Position_ID")

# Parse Position_ID format: "chr:start-end"
pos[, c("chr", "start_bp", "end_bp") := 
      tstrsplit(Position_ID, "[:-]", type.convert = TRUE, keep = 1:3)]

# Calculate locus length
pos[, length_bp := end_bp - start_bp + 1L]

# Sort by genomic position
setorder(pos, chr, start_bp)

# --- Format and Save Output ---
# Rename column for consistency
setnames(pos, "Position_ID", "pos_id")

# --- Save Outputs ---

if(!dir.exists("data/processed/joinLoci/relational_tables")){
  dir.create("data/processed/joinLoci/relational_tables", recursive = TRUE)
}

# Ensure the output directory exists
summary_dir <- dirname(opt$options$output_summary)
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
}

# Scans and thresholds are output in same directory
results_dir <- dirname(opt$options$output_scans)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

write.csv(merged_regions, opt$options$output_summary, row.names = FALSE)
write.csv(trait_loci, opt$options$output_relational_summary, row.names = FALSE)

saveRDS(scans, opt$options$output_scans)
saveRDS(thresholds, opt$options$output_thresholds)

# Save position reference table
fwrite(pos, opt$options$output_pos_summary)

# Print summary statistics
cat("Unique positions found:", nrow(pos), "\n")
cat("Chromosomes covered:", length(unique(pos$chr)), "\n")

