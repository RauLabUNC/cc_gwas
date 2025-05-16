# Plot allele effects
current_locus_info
pos_in_range_logical <- current_scan_data$pos$Mb * 10^6 >= plot_start_bp & current_scan_data$pos$Mb * 10^6 <= plot_end_bp
chr_logical <- current_scan_data$chr == current_locus_info$chr
marker_positions_bp <- current_scan_data$pos$Mb[pos_in_range_logical & chr_logical] * 10^6

### Fix this to subset the loci first by chromsome
markers <- data.frame(marker = current_scan_data$loci[pos_in_range_logical & chr_logical],
                      start  = marker_positions_bp) |> 
            arrange(start)

allele_effects_matrix <- current_scan_data$allele.effects[,markers$marker, drop = FALSE] # drop=FALSE to handle single marker case
allele_effects_transposed <- t(current_scan_data$allele.effects) |> as.data.frame()
allele_effects_transposed$marker <- rownames(allele_effects_transposed)


num_strains <- ncol(allele_effects_transposed) -1 
chromosome_name <- paste0("chr", current_locus_info$chr)

plot_data_list <- lapply(1:num_strains, function(i){
  curr_strain <- colnames(allele_effects_transposed)[[i]]
  temp_df <- allele_effects_transposed |> 
    dplyr::select(marker, founder_strains[i]) %>% 
    filter(marker %in% markers$marker)
  
  temp_df <- left_join(temp_df, markers, by = "marker")
  
  temp_df <- temp_df |> 
    mutate(chrom = chromosome_name) |> 
    arrange(start)
  colnames(temp_df)[2] <- "score"

  temp_df$end <- c(temp_df$start[2:nrow(temp_df)] - 1L, 
                   temp_df$start[nrow(temp_df)] + 1)   
  temp_df <- temp_df |> 
    dplyr::select(chrom, start, end, score)
  return(temp_df)
  }
)
names(plot_data_list) <- founder_strains



strain_colors <- rep(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), length.out = num_strains)

page_width_inches <- 6
page_height_inches <- 4 
plot_region_x <- 0.5
plot_region_y <- 0.5
plot_region_width <- 5
plot_region_height <- 2.5 
legend_y_start <- plot_region_y + plot_region_height + 0.2 
y_range <- c(-1,1)

plotgardener::pageCreate(width = page_width_inches, height = page_height_inches, default.units = "inches")

params_region <- plotgardener::pgParams(
  chrom = chromosome_name,
  chromstart = plot_start_bp,
  chromend = plot_end_bp,
  assembly = mm39
)

first_strain_name <- founder_strains[1]
for(i in 1:length(plot_data_list)){
  plotgardener::plotSignal(
    data = plot_data_list[[i]],
    params = params_region,
    range = y_range, 
    linecolor = strain_colors[i],
    fill = NA, 
    x = plot_region_x,
    y = plot_region_y,
    width = plot_region_width,
    height = plot_region_height,
    just = c("left", "top"),
    default.units = "inches",
    baseline = TRUE, 
    baseline.color = "grey"
  )
}
