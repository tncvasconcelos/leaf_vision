library(dplyr)
library(ggplot2)
# setwd("~/leaf_vision/")

merged_dataset <- read.csv("data/merged_dataset.csv")
merged_dataset$LMA <- merged_dataset$LMA*100

# Calculate mean LMA per 1 degree of latitude
mean_LMA_by_lat <- merged_dataset %>%
  filter(!is.na(LMA), !is.na(lat)) %>%  # Remove rows with missing values
  mutate(lat_bin = round(lat)) %>%  # Round latitude to the nearest degree
  group_by(lat_bin) %>%  # Group by latitude bin
  summarize(mean_LMA = mean(LMA, na.rm = TRUE), .groups = "drop")  # Calculate mean LMA

# Plot mean LMA as a line plot
mean_lma_ldg <- ggplot(mean_LMA_by_lat, aes(x = lat_bin, y = log(mean_LMA))) +
  geom_line(color = "forestgreen", size = 0.5) +  # Line plot
  theme_minimal() +
  labs(title = "Mean LMA by Latitude",
       x = "Latitude (degrees)", 
       y = "Mean LMA") +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  coord_flip()

# Count the number of observations per 1-degree latitude
counts_by_lat <- merged_dataset %>%
  filter(!is.na(lat)) %>%  # Remove rows with missing latitude values
  mutate(lat_bin = round(lat)) %>%  # Round latitude to the nearest degree
  group_by(lat_bin) %>%  # Group by latitude bin
  summarize(count = n(), .groups = "drop")  # Count the number of observations

# Plot number of observations as a line plot
count_ldg <- ggplot(counts_by_lat, aes(x = lat_bin, y = log(count))) +
  geom_line(color = "purple3", size = 0.5) +  # Line plot
  theme_minimal() +
  labs(title = "Number of Observations by Latitude",
       x = "Latitude (degrees)",
       y = "Number of Observations") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))+
  coord_flip()

pdf("FIGURES/plot_ldg.pdf" ,height=10,width=2)
grid.arrange(mean_lma_ldg, count_ldg, ncol=1, nrow = 2)
dev.off()



# # Plot LMA vs. Latitude with a LOESS line
# ggplot(merged_dataset, aes(x = lat, y = LMA)) +
#   geom_point(alpha = 0.5, color = "skyblue", fill="red") +  # Scatter points with transparency
#   geom_smooth(method = "loess", color = "red", se = TRUE) +  # LOESS line with confidence interval
#   theme_minimal() +
#   labs(title = "LOESS Smoothing of LMA by Latitude",
#        x = "Latitude",
#        y = "LMA") +
#   theme(panel.grid = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12))
