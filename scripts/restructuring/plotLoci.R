# This script analyzes significant QTL regions identified from different
# normalization methods (Box-Cox vs. Z-Score). It identifies potentially
# pleiotropic regions and compares the impact of the normalization methods
# by plotting LOD scores for significant loci alongside alternative conditions.

# Load Libraries
library(tidyverse)
library(plotgardener)
library(TxDb.Mmusculus.UCSC.mm39.knownGene) # For mm39
library(org.Mm.eg.db)                     

# Load Data (Assuming these are loaded in the environment as per the prompt)
final_sig_regions_df <- read.csv("results/sig_regions/all_significant_regions_summary.csv")
scan_data <- readRDS("results/sig_regions/scan_data.rds")
threshold_data <- readRDS("results/sig_regions/threshold_data.rds")


# Ensure position columns are numeric and filter invalid rows
sig_regions_df <- final_sig_regions_df %>%
  mutate(across(c(upper_pos_lod_drop, peak_pos, lower_pos_lod_drop, max_lod), as.numeric)) %>%
  filter(!is.na(upper_pos_lod_drop) & !is.na(lower_pos_lod_drop) & !is.na(chr) & !is.na(trait) & !is.na(norm) & !is.na(drug))

# Load genome assembly 
mm39 <- assembly(Genome = "mm39_GRCm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene", OrgDb = "org.Mm.eg.db")

# Define output directory
output_base_dir <- "results/sig_loci_plots_comparison_v3" # Consider updating version in name if desired
if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir, recursive = TRUE)
}


# Define function to get and format scan data (User Provided, slightly modified for safety)
get_scan <- function(n, t, d, current_chr_val) {
  scan.df <- NULL # Initialize as NULL
  scan_raw <- tryCatch(
    scan_data[[n]][[t]][[d]]
  )
  
    # Ensure marker names exist, create if necessary
    marker_names <- names(scan_raw$LOD)
    
    scan.df <- data.frame(marker = marker_names,
                          chr = as.character(scan_raw$chr), # Ensure character
                          pos = scan_raw$pos$Mb,
                          lod = scan_raw$LOD,
                          stringsAsFactors = FALSE) # Avoid factors
    
    # Convert position and filter
    scan.df <- scan.df %>%
      mutate(
        pos = as.numeric(pos) * 1e6 # Convert Mb to bp
      ) %>%
      filter(chr == current_chr_val & !is.na(pos) & !is.na(lod)) # Filter by current chromosome and remove NA pos/lod
  return(scan.df)
}


# Iterate through each significant region
for (i in 1:nrow(sig_regions_df)) {
  sig_row <- sig_regions_df[i, ]
  
  # Extract key info for the significant locus
  current_norm <- sig_row$norm
  current_trait <- sig_row$trait
  current_drug <- sig_row$drug
  current_chr <- as.character(sig_row$chr)
  # Add "chr" prefix if not present (adjust if your chromosome names differ)
  current_chr_str <- if (!startsWith(current_chr, "chr")) paste0("chr", current_chr) else current_chr
  
  peak_pos_bp <- sig_row$peak_pos * 1e6
  lower_bound_bp <- sig_row$lower_pos_lod_drop * 1e6
  upper_bound_bp <- sig_row$upper_pos_lod_drop * 1e6
  
  # Determine alternative conditions
  alt_drug <- if (current_drug == "Ctrl") "Iso" else "Ctrl"
  alt_norm <- if (current_norm == "boxcox") "zscore" else "boxcox"
  
  # Define plot region based on LOD drop interval
  chromstart <- floor(min(lower_bound_bp, upper_bound_bp)) # Use min/max for safety
  chromend <- ceiling(max(lower_bound_bp, upper_bound_bp))
  
  # Add padding around the interval for better visualization
  interval_width <- chromend - chromstart
  # Ensure interval_width is valid before calculating padding
  if (is.na(interval_width) || interval_width <= 0) interval_width <- 1 # Avoid issues with zero/negative width
  padding <- max(500000, interval_width * 0.2) # Add at least 500kb or 20% padding
  chromstart <- floor(max(0, chromstart - padding)) # Ensure start is not negative
  chromend <- ceiling(chromend + padding)
  
  
  # --- Retrieve Scan Data ---
  scan_original <- get_scan(current_norm, current_trait, current_drug, current_chr)
  scan_alt_drug <- get_scan(current_norm, current_trait, alt_drug, current_chr)
  scan_alt_norm <- get_scan(alt_norm, current_trait, current_drug, current_chr)
  scan_alt_norm_drug <- get_scan(alt_norm, current_trait, alt_drug, current_chr)
  
  # --- Retrieve Threshold Data ---
  get_thresh <- function(n, t, d) {
    thresh <- tryCatch(
      threshold_data[[n]][[t]][[d]],
      error = function(e) NULL # Return NULL on error
    )
    # Ensure threshold is a single numeric value
    if (!is.null(thresh) && (!is.numeric(thresh) || length(thresh) != 1)) {
      warning(paste("Threshold for", n, t, d, "is not a single numeric value. Setting to NULL."))
      thresh <- NULL
    }
    return(thresh)
  }
  
  thresh_original <- get_thresh(current_norm, current_trait, current_drug)
  thresh_alt_drug <- get_thresh(current_norm, current_trait, alt_drug)
  thresh_alt_norm <- get_thresh(alt_norm, current_trait, current_drug)
  thresh_alt_norm_drug <- get_thresh(alt_norm, current_trait, alt_drug)
  
  # Skip if original scan data is missing or has no rows after filtering
  if (is.null(scan_original) || nrow(scan_original) == 0) {
    warning(paste("Skipping row", i, ": Original scan data missing, empty, or invalid structure for",
                  current_norm, current_trait, current_drug, "on chr", current_chr))
    next # Skip to the next iteration
  }
  
  # --- Determine Y-axis limits ---
  max_lod_values <- c(0) # Start with 0 to ensure axis starts at zero
  if (!is.null(scan_original) && nrow(scan_original) > 0) max_lod_values <- c(max_lod_values, max(scan_original$lod, na.rm = TRUE))
  if (!is.null(scan_alt_drug) && nrow(scan_alt_drug) > 0) max_lod_values <- c(max_lod_values, max(scan_alt_drug$lod, na.rm = TRUE))
  if (!is.null(scan_alt_norm) && nrow(scan_alt_norm) > 0) max_lod_values <- c(max_lod_values, max(scan_alt_norm$lod, na.rm = TRUE))
  if (!is.null(scan_alt_norm_drug) && nrow(scan_alt_norm_drug) > 0) max_lod_values <- c(max_lod_values, max(scan_alt_norm_drug$lod, na.rm = TRUE))
  
  # Include non-NULL, numeric thresholds in the max calculation
  threshold_values <- na.omit(c(thresh_original, thresh_alt_drug, thresh_alt_norm, thresh_alt_norm_drug))
  # Filter again to ensure they are numeric (get_thresh should handle this, but belt-and-suspenders)
  threshold_values <- threshold_values[sapply(threshold_values, is.numeric)]
  if (length(threshold_values) > 0) {
    max_lod_values <- c(max_lod_values, threshold_values)
  }
  
  max_observed <- max(max_lod_values, na.rm = TRUE)
  # Handle case where max_observed might still be -Inf if all data is NA/NULL
  if (!is.finite(max_observed)) max_observed <- 0
  ylim_max <- ceiling(max_observed) # Add 10% buffer and ensure at least 1 unit above max
  ylim_max <- max(ylim_max, 5) # Ensure minimum ylim max of 5 if max_observed is low
  plot_ylim <- c(0, ylim_max) # Final y-axis limits for the plot
  
  # --- Create Plot ---
  output_filename <- file.path(output_base_dir,
                               paste0("locus_", current_norm, "_", current_trait, "_",
                                      current_drug, "_chr", current_chr, "_peak",
                                      round(sig_row$peak_pos, 2), "Mb.png"))
  
# Define plot dimensions and positions
  plot_x <- 3.25
  plot_y <-  0.5# Top y-coordinate of the signal plot area
  plot_width <- 6
  plot_height <- 2.2
  gene_height <- 2 # Adjusted height for genes
  gene_y_offset <- 0.25 # Space between signal plot and genes
  label_y_offset <- 0.1 # Space between genes and genome label
  
  
  # Define common parameters for the genomic region and assembly
  # Select the correct assembly object based on your data
  current_assembly <- mm39 # Or mm39 if using mouse data
  params_genome <- pgParams(
    assembly = current_assembly,
    chrom = current_chr_str,
    chromstart = chromstart,
    chromend = chromend
  )
  

png(file = output_filename, width = 6.5, height = 7.5 , units = "in", res = 300)
  
  # Create a plotgardener page
  pageCreate(width = 6.5, height = 7, default.units = "inches", showGuides = F)
  
  # 1. Plot each LOD curve and its threshold line sequentially
  plot_lod_curve <- function(scan_df, drug, norm, threshold_val, type = "", is_first_plot = FALSE) {
    #scan_df <- scan_original
    #drug <- current_drug
    #norm <- current_norm
    #type = ""
    #threshold_val <- thresh_original
    scan_obj <- scan_df |> 
      mutate(chrom = paste0("chr", chr),
             p = 10^(-lod)) |> 
      dplyr::select(chrom, pos, p)  
    y <- ifelse(type == "alt_drug", plot_y+plot_height+0.2, plot_y)
      # Plot LOD score data using plotSignal
      manhattanPlot <- plotManhattan(
        data = scan_obj,
        params = params_genome,
        trans = "-log10", sigVal = 10^(-threshold_val),
        range = plot_ylim,
        x = 3.25, y = y, width = 6, height = plot_height,
        just = c("center", "top"), 
        xfield = "pos", yfield = "p", fill = colors[paste0(norm, ".nonSig")], sigCol =  colors[paste0(norm, ".Sig")],
        sigLine = T, baseline = T,
        yscale_reverse = F
      )
    return(manhattanPlot) # Return the plot object if it was the first, otherwise NULL
  }
  # Define color pal
  colors <- c("boxcox" =  c(Sig = "#1f78b4" , nonSig = "#a6cee3"), "zscore" = c(Sig = "#ff7f00" , nonSig = "#fdbf6f"))
  # Plot the curves - plot the one guaranteed to exist (original) first
  # to establish the viewport and get a reference object.
  main.plot <- plot_lod_curve(scan_original, drug = current_drug,
               norm = current_norm, threshold_val = thresh_original, is_first_plot = TRUE)
  annoYaxis(plot = main.plot,
    at = pretty(plot_ylim), # Use pretty() for nice tick marks
    axisLine = TRUE, fontsize = 8, main = FALSE # Ensure main=FALSE
  )
  # Plot other curves (they will overlay the first one)
  plot_lod_curve(scan_alt_norm, drug = current_drug,
                 norm = alt_norm, threshold_val = thresh_alt_norm)
  sub.plot <- plot_lod_curve(scan_alt_norm_drug, drug = alt_drug,
                 norm = alt_norm, threshold_val = thresh_alt_norm_drug, type = "alt_drug")
  plot_lod_curve(scan_alt_drug, drug = alt_drug,
                 norm = current_norm, threshold_val = thresh_alt_drug, type = "alt_drug")
  
  annoYaxis(plot = sub.plot,
            at = pretty(plot_ylim), # Use pretty() for nice tick marks
            axisLine = TRUE, fontsize = 8, main = FALSE # Ensure main=FALSE
  )
  # 2. Annotate the Y axis onto the reference plot viewport

  # Add Y-axis label text separately using plotText
  plotText(label = paste0("LOD — ", current_drug),
           x = 2*plot_x + 0.1, y = plot_y + plot_height/2, # Position left of the plot area
           rot = 270, # Rotate text
           fontsize = 10, just = "center",
           default.units = "inches")
  
  # Add Y-axis label text separately using plotText
  plotText(label = paste0("LOD — ", alt_drug),
           x = 2*plot_x + 0.1, y = plot_height/2 + 2.4 + plot_y, # Position left of the plot area
           rot = 270, # Rotate text
           fontsize = 10, just = "center",
           default.units = "inches")
  
  # 3. Highlight original significant region interval using annoHighlight
  # Reference the established plot viewport
  annoHighlight(
    plot = main.plot, # Reference the plot object
    chrom = current_chr_str, # Already in params_genome, but explicit here
    chromstart = floor(min(lower_bound_bp, upper_bound_bp)), # Use the actual interval bounds
    chromend = ceiling(max(lower_bound_bp, upper_bound_bp)),
    fill = "#fb9a99", # Semi-transparent gold
    y = plot_y, height = plot_height, # Match the signal plot viewport y and height
    just = c("left", "top"), # Anchor point
    default.units = "inches",
    alpha = 0.2,
    params = params_genome # Pass genomic parameters
  )
  
  # 5. Plot Genes below the signal plot
  gene_y_pos <- 4.5 + plot_y # Calculate y position
  # Assign the plotGenes output to potentially use for annoGenomeLabel
  genes_plot <- plotGenes(
    params = params_genome, # Pass genomic parameters
    x = plot_x, y = gene_y_pos, width = 6, height = gene_height, # Define viewport
    stroke = 1, fontsize = 10, # Styling
    just = c("center", "top"), # Anchor point
    default.units = "inches"
  )
  
  # 6. Plot Genome Label below genes
  label_y_pos <- gene_y_pos + gene_height  # Calculate y position
  annoGenomeLabel(
    plot = genes_plot, # Use the genes plot object as reference
    params = params_genome, # Pass genomic parameters (redundant if using plot obj, but safe)
    x = plot_x, y = label_y_pos, # Position label explicitly below genes
    scale = "Mb", fontsize = 8, # Styling
    just = c("center", "top"), # Anchor point
    default.units = "inches"
  )
  
  # 7. Add Title using plotText
  title_text <- paste("QTL:", current_trait, " (", current_norm, "/", current_drug, ") - ",
                      current_chr_str, " Peak:", round(sig_row$peak_pos, 2), "Mb")
  plotText(
    label = title_text,
    x = 3.25, y = 0, # Position near top-center of the page
    just = "center", fontface = "bold", fontsize = 12,
    default.units = "inches"
  )
  
  # 8. Add Legend using plotLegend
  # Only plot legend if there are items to include
    plotLegend(
      legend = c("Box-Cox", "Z-score"),
      fill = colors[c(1,3)], # Use fill for the color squares
      border = FALSE, # No border around legend items
      x = 3.25, # Position to the right of the signal plot
      y = plot_height+plot_y-0.2, # Align top with signal plot
      width = 1.0, height = 0.75, # Adjust size as needed
      just = c("center", "top"), # Anchor point
      fontsize = 10,
      default.units = "inches",
      orientation = "h"
    )
  
  
  # Close the PNG device
  dev.off()
}
print("Finished plotting significant loci comparisons.")
