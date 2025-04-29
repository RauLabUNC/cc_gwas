library(miqtl)


# Load the data
boxcox_hr <- readRDS("data/processed/ropscan/boxcox_individual_HR.0_Ctrl.rds")
zscore_hr <- readRDS("data/processed/ropscan/zscore_individual_HR.0_Ctrl.rds")
boxcox_bw <- readRDS("data/processed/ropscan/boxcox_individual_BW.day.0_Ctrl.rds")
zscore_bw <- readRDS("data/processed/ropscan/zscore_individual_BW.day.0_Ctrl.rds")


# put the p-values into tables

# Extract p-values and create dataframes for each dataset
boxcox_hr_df <- data.frame(
  SNP = names(boxcox_hr$p.value),
  p_value = boxcox_hr$p.value,
  trait = "HR.0",
  norm_method = "boxcox"
)

zscore_hr_df <- data.frame(
  SNP = names(zscore_hr$p.value),
  p_value = zscore_hr$p.value,
  trait = "HR.0",
  norm_method = "zscore"
)

boxcox_bw_df <- data.frame(
  SNP = names(boxcox_bw$p.value),
  p_value = boxcox_bw$p.value,
  trait = "BW.day.0",
  norm_method = "boxcox"
)

zscore_bw_df <- data.frame(
  SNP = names(zscore_bw$p.value),
  p_value = zscore_bw$p.value,
  trait = "BW.day.0",
  norm_method = "zscore"
)

# Combine all dataframes into one
combined_df <- bind_rows(boxcox_hr_df, zscore_hr_df, boxcox_bw_df, zscore_bw_df) |> 
  pivot_wider(values_from = p_value, names_from = norm_method)


scatterPlot <- combined_df |> 
  ggplot(aes(x = -log10(zscore), y = -log10(boxcox), color = trait)) +
    geom_abline(slope = 1, color = "black") +
    geom_point(alpha = 0.3) +
    facet_grid(~trait) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal()

ggsave("results/boxCox/hr_bw_zscoreVSboxcox.png", plot = scatterPlot,
       height = 4, width = 6, dpi = 200)

png(filename ="results/boxCox/hr_scan_zscoreVSboxcox.png", res = 300, height = 6, width = 9, units = "in")
genome.plotter.whole(scan.list=list(zscore_individuals = zscore_hr,
                                    boxcox_individuals = boxcox_hr),
                    main = "miQTL for HR.0, individuals mapped")

dev.off()

png(filename ="results/boxCox/bw_scan_zscoreVSboxcox.png", res = 300, height = 6, width = 9, units = "in")
genome.plotter.whole(scan.list=list(zscore_mean = zscore_bw,
                                    boxcox_mean = boxcox_bw),
                     main = "miQTL for BW.0, individuals mapped")

dev.off()

