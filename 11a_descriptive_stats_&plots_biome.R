

# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
library(gridExtra)
library(data.table)
library(ggplot2)

#---------------------------------------
merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset$LMA <- merged_dataset$LMA * 100
merged_dataset$LMA <- log(merged_dataset$LMA)

library(ggplot2)
library(ggridges)
library(viridis)

# Remove NAs from LMA and reorder biome by median LMA
merged_dataset <- na.omit(merged_dataset[ , c("LMA", "super_biome")])

# Create the ridgeline plot with median lines on each ridge
ggplot(merged_dataset, aes(x = LMA, y = super_biome, fill = super_biome)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, 
                      scale = 1.5, rel_min_height = 0.01, alpha = 0.5) +  # Add median line within each ridge
  scale_fill_viridis_d() +  # Viridis color palette for discrete values
  theme_bw() +  # Black-and-white theme
  labs(x = "LMA", y = "Biome", fill = "Biome")  # Axis and legend labels
