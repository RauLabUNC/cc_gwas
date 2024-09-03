# This script is meant to clean up the CC phenotype data file from August 2024

# Load libs 
library(tidyverse)

# Load data
data <- read.csv("data/raw/phenotypes/full_cc_panel_data_08_06_24.csv")

colnames(data)[which(colnames(data) == "Iso.Ctrl")] <- "Drug"

# Clean the strain names
data$Strain <- lapply(str_split(data$Strain, "[.]"), "[[", 1) |>  unlist()

data <- data |> 
  mutate(
    Strain_Clean = case_when(
      Strain %in% c(1:9) ~ paste0("CC00", Strain),
      Strain %in% c(10:99) ~ paste0("CC0", Strain),
      .default = Strain
    )
  )

data$Strain_Clean <- lapply(str_split(data$Strain_Clean, "/"), "[[", 1) |>  unlist()


# Manually match to the genotype cache files
data <- data |> 
  mutate(Strain_Clean = case_when(
    Strain_Clean == "A" ~ "AJ",
    Strain_Clean == "C57B" ~ "B6",
    Strain_Clean == "Cast" ~ "CAST",
    .default = Strain_Clean
  ))

# Check the names used in the genotypes 
load("data/raw/genomes/haplotype_cache_cc_083024/additive/chr1/subjects.RData")

# now, clean the sex and treatment values

data <- data |> 
  mutate(Drug_Clean = case_when(
    Drug == "Ctrl" ~ 0,
    Drug == "Iso" ~ 1),
         Sex_Clean = case_when(
      Sex == "M" ~ 0,
      Sex == "F" ~ 1
    ))

# Average them by strain, sex, and treatment
# Group phenotypes by strain, drug, and sex, then calculate mean of phenotype columns
grouped_phenotypes <- data %>%
  group_by(Strain_Clean, Drug_Clean, Sex_Clean) %>%
  summarize(across(BW.day.0:delta.EF, ~mean(as.numeric(.), na.rm = TRUE)))

# Ensure the directory exists
output_dir <- "data/processed/phenotypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(grouped_phenotypes, file.path(output_dir, "mean_cc_panel_08_06_24.csv"), row.names = F)