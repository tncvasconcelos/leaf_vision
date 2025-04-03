# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
# setwd("~/leaf_vision/")
library(gridExtra)
library(data.table)
library(ggplot2)

# gt petiole comparison

manual_measurements <- as.data.frame(fread("data/GT_comparison.csv"))
load("results/data_subset_for_plots.Rsave")
manual_measurements <- merge(manual_measurements,
                             dat,
                             by.x = "filepart1", by.y="sp")

manual_measurements$lma <- exp(manual_measurements$lma)*100
#model <- lm(manual_measurements$width_pixels~manual_measurements$pixel_distance)
manual_measurements$abs_error <- abs(manual_measurements$width_pixels-manual_measurements$pixel_distance)
model <- lm(manual_measurements$lma~manual_measurements$abs_error)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

comparison_scatter_plot <- ggplot(manual_measurements, aes(x = abs_error, y = lma)) +
  geom_point(aes(color = lma), size = 3, alpha = 0.7) +
  scale_color_viridis_c(option = "C", end = 0.85) +
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(manual_measurements$abs_error) , y = max(manual_measurements$lma) * 0.9, 
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    y = "LMA",
    x = "absolute error petiole width",
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

pdf("plots/LMAvsError.pdf")
comparison_scatter_plot
dev.off()
