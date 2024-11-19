#LMA 
library(scales)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

setwd("~/leaf_vision/")

merged_dataset <- read.csv("data/merged_dataset.csv")
merged_dataset$LMA <- merged_dataset$LMA*100

lma_results <- aggregate(merged_dataset$LMA, list(merged_dataset$filename), 
                         FUN = function(x) c(mean(x)))
merged_dataset <- subset(merged_dataset, !duplicated(merged_dataset$filename))
merged_dataset <- merge(merged_dataset, lma_results, by.x="filename", by.y="Group.1")
merged_dataset$LMA <- merged_dataset$x

#-------------------
# MAPS WITH GLOBAL DISTRIBUTION
world <- ne_countries(scale = "medium", returnclass = "sf")

#-------------------
merged_dataset$grid_lat <- floor(merged_dataset$lat / 2.5) * 2.5
merged_dataset$grid_lon <- floor(merged_dataset$lon / 2.5) * 2.5

# Aggregate data
grid_data <- merged_dataset %>%
  group_by(grid_lat, grid_lon) %>%
  summarize(mean_LMA = mean(LMA, na.rm = TRUE), .groups = "drop")

mean_lma_map <- ggplot(data = world) +
  geom_sf(fill = "lightgrey", color = "white") +  # World map base
  geom_tile(data = grid_data, 
            aes(x = grid_lon, y = grid_lat, fill = log(mean_LMA)), 
            width = 2.5, height = 2.5, alpha = 0.8) +  # 5x5 grid cells
  scale_fill_viridis_c(option = "viridis", name = "Mean LMA") +  # Color scale
  theme_minimal() +
  labs(title = "Mean LMA by 2.5x2.5 Degree Grid",
       x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))

#-------------------
merged_dataset$grid_lat <- floor(merged_dataset$lat / 2.5) * 2.5
merged_dataset$grid_lon <- floor(merged_dataset$lon / 2.5) * 2.5

grid_density <- merged_dataset %>%
  group_by(grid_lat, grid_lon) %>%
  summarize(count = n(), .groups = "drop")

# Plot the heatmap
n_points_map <- ggplot(data = world) +
  geom_sf(fill = "lightgrey", color = "white") +  # World map base
  geom_tile(data = grid_density, 
            aes(x = grid_lon, y = grid_lat, fill = log(count)), 
            width = 2.5, height = 2.5, alpha=0.8) +  # 5x5 grid cells
  scale_fill_viridis_c(option = "plasma", name = "Point Count") +  # Density color scale
  theme_minimal() +
  labs(title = "Point Density Heatmap by 2.5x2.5 Degree Grid",
       x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))

pdf("FIGURES/plot_maps.pdf" ,height=10,width=10)
grid.arrange(mean_lma_map, n_points_map, ncol=1, nrow = 2)
dev.off()

#------------------
# VIOLIN PLOTS PER BIOME

# Filter dataset and remove biomes with fewer than 100 observations
filtered_dataset <- merged_dataset %>%
  filter(!is.na(LMA), !is.na(biome)) %>%
  group_by(biome) %>%
  filter(n() >= 100) %>%  # Only keep biomes with at least 100 observations
  ungroup()

# Reorder biomes by median LMA
filtered_dataset <- filtered_dataset %>%
  mutate(biome = reorder(biome, LMA, FUN = median, decreasing = TRUE))

summary_data <- filtered_dataset %>%
  group_by(biome) %>%
  summarize(
    count = n(),
    median_LMA = median(LMA, na.rm = TRUE),
    .groups = "drop"
  )

violin_biomes <- ggplot(data = filtered_dataset, aes(x = biome, y = LMA)) +
  geom_violin(fill = "skyblue", color = "darkblue", alpha = 0.7) +  # Violin plot
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violins
  geom_text(data = summary_data, 
            aes(x = biome, y = 500, label = count), 
            color = "black", vjust = -0.5, size = 3.5) +  # Add sample sizes above violins
  geom_text(data = summary_data, 
            aes(x = biome, y = median_LMA, label = round(median_LMA, 1)), 
            color = "red", vjust = -0.5, size = 3.5) +  # Add medians above violins
  coord_cartesian() +  # Limit the flipped axis (x-axis now represents LMA)
  theme_minimal() +
  labs(title = "",
       x = "", y = "LMA") +
  coord_flip(ylim = c(0, 500))
  

pdf("FIGURES/violin_plots_biomes_100_or_more.pdf", height=6, width=8)
violin_biomes
dev.off()

#---------------------------
