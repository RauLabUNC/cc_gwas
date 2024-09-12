# This script is meant to take the 'cleaned' phenotypes from 09/05/2024 and..
# .. format them for miQTL scans. Outliers have been noted and some of them..
# .. have been removed manually by CR. Todd needs to check the rest.

# Load libraries
library(xlsx)
library(tidyverse)

# Load data
phenos <- read.xlsx("data/raw/phenotypes/phenotypes_outliers_fibrosis_091224.xlsx", sheetIndex = 1)[,1:63]


# Create a function to ID outliers and convert them to NAs, for each phenotype
remove_outliers <- function(data, variables) {
  for (var in variables) {
    # Calculate the lower and upper bounds for outliers
    q1 <- quantile(data[[var]], 0.25, na.rm = T)
    q3 <- quantile(data[[var]], 0.75, na.rm = T)
    iqr <- q3 - q1
    lower_bound <- q1 - 3 * iqr
    upper_bound <- q3 + 3 * iqr
    print(var)
    # Remove outliers
    # Needed to not be NA when assigning values
    data[((!is.na(data[[var]]) & data[[var]] < lower_bound ) | 
          (!is.na(data[[var]]) & data[[var]] > upper_bound)), var] <- NA
  }
  
  return(data)
}

phenos <- phenos[,c(1:61,63,62)]
# List variables 
vars <- colnames(phenos)[10:62]
# Remove outliers
data.no.outliers <- remove_outliers(phenos, vars)

# Factorize variables + take means
data.no.outliers <- data.no.outliers |> 
  mutate(Drug_Clean = factor(case_when(
    Drug == "Ctrl" ~ 0,
    Drug == "Iso" ~ 1)),
    Sex_Clean = factor(case_when(
      Sex == "M" ~ 0,
      Sex == "F" ~ 1)),
    pheno.id = paste(Strain_Clean, Sex_Clean, Drug_Clean, sep = "_")) 

# Ensure the directory exists
output_dir <- "data/processed/phenotypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(data.no.outliers, file.path(output_dir, "no_outliers_cc_panel_08_06_24.csv"), row.names = F)
