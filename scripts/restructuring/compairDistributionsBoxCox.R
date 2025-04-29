# Load required libraries
library(tidyverse)

boxcox  <- read.csv("data/processed/phenotypes/boxCoxTest/boxcox_individual.csv")
zscores <- read.csv("data/processed/phenotypes/boxCoxTest/zscore_individual.csv")

boxcox <- boxcox |> dplyr::select(HR.0, BW.day.0, Drug) |> mutate(transform = "Box-Cox")

zscores <- zscores |> dplyr::select(HR.0, BW.day.0, Drug) |> mutate(transform = "Z-Score")

merged <- rbind(boxcox, zscores) |> pivot_longer(cols = c("HR.0", "BW.day.0"), names_to = "Trait")


merged |> 
  ggplot(aes(x = value, fill = transform)) +
  geom_density(alpha = 0.6) +
  facet_grid(transform~Trait, scales = "free_x")


boxcox |> 
  ggplot(aes(x = HR.0, fill = transform)) +
  geom_density(alpha = 0.6) 
zscores |> 
  ggplot(aes(x = HR.0, fill = transform)) +
  geom_density(alpha = 0.6) 


cor(boxcox$HR.0, zscores$HR.0, use = "complete")

cor(boxcox$HR.0, zscores$HR.0, use = "complete")

colnames(boxcox)
