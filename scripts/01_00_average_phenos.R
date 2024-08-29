# This script is meant to convert the phenotype measures from the CC HF...
# project from individuals to group-level means

# Load libraries
library(tidyverse)

# Load data
phenotypes <- read.csv("data/raw/phenotypes/full_cc_panel_data_04_16_24.csv")

# Group phenotypes by strain, drug, and sex, then calculate mean of phenotype columns
grouped_phenotypes <- phenotypes %>%
  group_by(Strain_Clean, Drug_Clean, Sex_Clean) %>%
  summarize(across(BW.day.0:delta.EF, ~mean(as.numeric(.), na.rm = TRUE)))

# Ensure the directory exists
output_dir <- "data/processed/phenotypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(grouped_phenotypes, file.path(output_dir, "mean_cc_panel_04_16_24.csv"), row.names = F)
