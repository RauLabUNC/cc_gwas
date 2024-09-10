# This script is meant to identify outlier values in the phenotype data

# Load data
# Load libs 
library(tidyverse)
library(data.table)
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

# now, clean the sex and treatment values


detectOutliers <- function(var){
    # Calculate the lower and upper bounds for outliers
    q1 <- quantile(data[[var]], 0.25, na.rm = T)
    q3 <- quantile(data[[var]], 0.75, na.rm = T)
    iqr <- q3 - q1
    lower_bound <- q1 - 2 * iqr
    upper_bound <- q3 + 2 * iqr
    outs <- data[[var]] > upper_bound | data[[var]] < lower_bound
    return(outs)}


phenos <- colnames(data)[10:61]
outliers <- lapply(phenos, function(x){detectOutliers(x) |> as.list()})

df <- outliers |> rbindlist() |> t() |> as.data.frame()
colnames(df) <- phenos

df[is.na(df)] <- FALSE

df <- cbind(data[,1:9], df)
df <- cbind(df, data$Strain_Clean)

library(xlsx)

# Create a new workbook and sheet
sheetname <- "mysheet"
write.xlsx(data, "data.xlsx", sheetName = sheetname, row.names = FALSE)
file <- "data.xlsx"

# Load the workbook
wb <- loadWorkbook(file)

# Create fill objects and cell styles
fo1 <- Fill(foregroundColor = "yellow", backgroundColor = "yellow", pattern = "SOLID_FOREGROUND")
cs1 <- CellStyle(wb, fill = fo1)

# Get the sheet and rows
sheets <- getSheets(wb)
sheet <- sheets[[sheetname]]
rows <- getRows(sheet, rowIndex = 2:(nrow(data) + 1))  # Adjust for header row

# Get the cells
cols <- ncol(data)
cells <- getCells(rows, colIndex = 1:cols)

# Extract the cell values
values <- lapply(cells, getCellValue)

# Find cells that need to be highlighted
highlightyellow <- NULL
for (i in seq_along(values)) {
  row <- (i - 1) %/% cols + 1
  col <- (i - 1) %% cols + 1
  if (df[row, col] == TRUE) {
    highlightyellow <- c(highlightyellow, names(values)[i])
  }
}

# Apply the formatting
lapply(names(cells[highlightyellow]), function(ii) setCellStyle(cells[[ii]], cs1))

# Save the workbook
saveWorkbook(wb, file)