# This script is meant to take the 'cleaned' phenotypes from 09/05/2024 and..
# .. format them for miQTL scans. Outliers have been noted and some of them..
# .. have been removed manually by CR. Todd needs to check the rest.

# Load libraries
library(xlsx)
library(tidyverse)

# Load data
phenos <- read.xlsx("data/raw/phenotypes/cc_phenotypes_outliers.xlsx", sheetIndex = 1)


# Create a function to ID outliers and convert them to NAs, for each phenotype
remove_outliers <- function(data, variables) {
  for (var in variables) {
    # Calculate the lower and upper bounds for outliers
    q1 <- quantile(data[[var]], 0.25, na.rm = T)
    q3 <- quantile(data[[var]], 0.75, na.rm = T)
    iqr <- q3 - q1
    lower_bound <- q1 - 2 * iqr
    upper_bound <- q3 + 2 * iqr
    print(var)
    # Remove outliers
    # Needed to not be NA when assigning values
    data[((!is.na(data[[var]]) & data[[var]] < lower_bound ) | 
          (!is.na(data[[var]]) & data[[var]] > upper_bound)), var] <- NA
  }
  
  return(data)
}

# List variables 
vars <- colnames(phenos)[10:61]
# Remove outliers
data.no.outliers <- remove_outliers(phenos, vars)

# Factorize variables + take means
data.no.outliers <- data.no.outliers |> 
  mutate(Drug_Clean = case_when(
    Drug == "Ctrl" ~ 0,
    Drug == "Iso" ~ 1),
    Sex_Clean = case_when(
      Sex == "M" ~ 0,
      Sex == "F" ~ 1
    ))

# Average them by strain, sex, and treatment
# Group phenotype by strain, drug, and sex, then calculate mean of phenotype columns
# Then convert them to z scores (subtract mean, divid by StdDev) to get to normal dist
grouped_phenotypes <- data.no.outliers |> 
  group_by(Strain_Clean, Drug_Clean, Sex_Clean) |> 
  summarize(across(BW.day.0:delta.EF, ~mean(as.numeric(.), na.rm = TRUE))) |>
  ungroup() |> 
  mutate(across(BW.day.0:delta.EF, ~scale(.)[, 1]))

# Ensure the directory exists
output_dir <- "data/processed/phenotypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(grouped_phenotypes, file.path(output_dir, "no_outliers_cc_panel_08_06_24.csv"), row.names = F)
