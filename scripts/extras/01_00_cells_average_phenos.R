# This script is meant to convert the phenotype measures from the CC HF...
# project from individuals to group-level means

# Load libraries
library(tidyverse)
library(xlsx)
library(compositions)

# Load data
phenotypes <- read.xlsx("data/raw/phenotypes/outliersTraits_03242025.xlsx", sheetIndex = 1)
#phenotypes <- read.csv("data/raw/phenotypes/outliersTraits_03312025.csv")
info <- read.csv("data/raw/phenotypes/sampleInfoBasic.csv")
props <- read.csv("data/processed/phenotypes/cell_comps_04022025.csv")

# Convert proportions to CLRs
props_wide <- props |> 
  pivot_wider(names_from = CellType, values_from = Prop, id_cols = CommonID)

clr_props <- sapply(1:nrow(props_wide), function(x){clr(props_wide[x,-1])}) |> 
  t() |> 
  as.data.frame()

clr_props$CommonID <- props_wide$CommonID

# Add sampleID from info to props
info$CommonID <- paste(info$PlateNumber, info$AlitheaBarcode, sep = "_")

info_sub <- info |> filter(PlateNumber %in% c(1,2,5) & Drug == "Iso")

props_join <- inner_join(info_sub, clr_props)

#props_join <- pivot_wider(props_join, names_from = CellType, values_from = Prop)
props_join <- props_join |>  mutate(Endothelial_Cells = `Endothelial Cells`) |> select(-`Endothelial Cells`)

all_join <- left_join(props_join, phenotypes)

# Factorize variables + take means
phenotypes_factor <- all_join |> 
  mutate(Drug_Binary = factor(case_when(
    Drug == "Ctrl" ~ 0,
    Drug == "Iso" ~ 1)),
    Sex_Binary = factor(case_when(
      Sex == "M" ~ 0,
      Sex == "F" ~ 1))) 

# Group phenotypes by strain, drug, and sex, then calculate mean of phenotype columns
scaled_phenotypes <- phenotypes_factor |> 
  mutate(gwas_temp_id = paste0(Strain, Drug_Binary, Sex_Binary)) |> 
  #group_by(gwas_temp_id, Strain, Drug_Binary, Sex_Binary) |> 
  #summarize(across(BW.day.0:CSA, ~mean(as.numeric(.), na.rm = TRUE))) |>
  group_by(Drug_Binary) |> # Center values to zero and scale within each treatment
  mutate(across(BW.day.0:CSA, ~ (. - mean(., na.rm = T))/sd(., na.rm = T))) |> 
  as.data.frame()


# Ensure the directory exists
output_dir <- "data/processed/phenotypes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.csv(scaled_phenotypes, file.path(output_dir, "meanCenterScaledByTreat_04022025.csv"), row.names = F)
