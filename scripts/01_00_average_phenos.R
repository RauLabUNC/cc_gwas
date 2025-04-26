# This script is meant to convert the phenotype measures from the CC HF...
# project from individuals to group-level means

# Load libraries
library(tidyverse)
library(xlsx)

# Load data
phenotypes <- read.xlsx("data/raw/phenotypes/outliersTraits_03242025.xlsx", sheetIndex = 1)


# Factorize variables + take means
phenotypes_factor <- phenotypes |> 
  mutate(Drug_Binary = factor(case_when(
    Drug == "Ctrl" ~ 0,
    Drug == "Iso" ~ 1)),
    Sex_Binary = factor(case_when(
      Sex == "M" ~ 0,
      Sex == "F" ~ 1))) 

# Group phenotypes by strain, drug, and sex, then calculate mean of phenotype columns
scaled_phenotypes <- phenotypes_factor |> 
  mutate(gwas_temp_id = paste0(Strain, Drug_Binary, Sex_Binary)) |> 
  group_by(gwas_temp_id, Strain, Drug_Binary, Sex_Binary) |> 
  summarize(across(BW.day.0:CSA, ~mean(as.numeric(.), na.rm = TRUE))) |>
  group_by(Drug_Binary) |> # Center values to zero and scale within each treatment
  mutate(across(BW.day.0:CSA, ~ (. - mean(., na.rm = T))/sd(., na.rm = T))) |> 
  as.data.frame()


# Ensure the directory exists
output_dir <- "data/processed/phenotypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(scaled_phenotypes, file.path(output_dir, "meanCenterScaledByTreat_03242025.csv"), row.names = F)
