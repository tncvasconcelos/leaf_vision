# rm(list=ls())
# setwd("~/leaf_computer_vision")
library(gridExtra)
library(data.table)
library(ggplot2)

#---------------------------------------
merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$area))
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$petiole_width))

# how many specimens:
length(unique(merged_dataset$genus_species))

# how many species:
length(unique(merged_dataset$filename))

# how many families:
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- subset(names_sample, names_sample$taxon_name %in% unique(merged_dataset$genus_species))
length(unique(names_sample$family))
write.csv(names_sample, file="supporting_datasets/species_list_w_families.csv")

#---------------------------------------
# Some summary stats and plots:
min(merged_dataset$area)
merged_dataset$filename[which.min(merged_dataset$area)]
max(merged_dataset$area)
merged_dataset$filename[which.max(merged_dataset$area)]

# pdf("results/leaf_area_dist.pdf")
# hist(log(merged_dataset$area),breaks=100, xlab="log(leaf area cm^2)", main="Leaf area distribution")
# dev.off()

min(merged_dataset$petiole_width)
merged_dataset$filename[which.min(merged_dataset$petiole_width)]
max(merged_dataset$petiole_width)
merged_dataset$filename[which.max(merged_dataset$petiole_width)]

# pdf("results/petiole_width_dist.pdf")
# hist(log(merged_dataset$petiole_width), breaks=100, xlab="log(petiole width cm)", main="Petiole width distribution")
# dev.off()

#pdf("results/cor_leaf_area_petiole_width.pdf")
model <- lm(log(merged_dataset$area)~log(merged_dataset$petiole_width))
summary(model)
plot(log(merged_dataset$area)~log(merged_dataset$petiole_width), 
     xlab="log(petiole width cm)", ylab="log(leaf area cm^2)")

abline(model, col="red")
#dev.off()

#-------------------------------
# Histogram for PW
# Calculate the median
merged_dataset$petiole_width <- merged_dataset$petiole_width * 10
median_value <- median(merged_dataset$petiole_width)

petiole_width_hist <- ggplot(merged_dataset, aes(x = petiole_width)) +
  geom_histogram(binwidth = 0.1, fill = "orange", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = median_value), color = "black", linetype = "dashed", size = 1,
             show.legend = TRUE) +
  labs(title = "",
       x = "Petiole width (mm)", 
       y = "Frequency") +
  annotate("text", x = median_value + 0.9, y = 1500, label = paste("Median:", round(median_value, 1)),
           vjust = -0.5, color = "black", size = 5) +
  # Set x-axis limits from 0 to 500
  xlim(0, 5) +
  ylim(0, 4500) +
  theme_bw(base_size = 15) 

#-------------------------------
# Histogram for LMA
# Calculate the median
merged_dataset$LMA <- merged_dataset$LMA*100
median_value <- median(merged_dataset$LMA)

LMA_hist <- ggplot(merged_dataset, aes(x = LMA)) +
  geom_histogram(binwidth = 10, fill = "orange", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = median_value), color = "black", linetype = "dashed", size = 1,
             show.legend = TRUE) +
  labs(title = "",
       x = expression("LMA " (g/m^2)), 
       y = "Frequency") +
  annotate("text", x = median_value + 60, y = 1500, label = paste("Median:", round(median_value, 1)),
           vjust = -0.5, color = "black", size = 5) +
  # Set x-axis limits from 0 to 500
  xlim(0, 500) +
  theme_bw(base_size = 15) 

#-------------------------------
# Histogram for area
# Calculate the median
merged_dataset$area <- merged_dataset$area * 0.0001
median_value <- median(merged_dataset$area)

LA_hist <- ggplot(merged_dataset, aes(x = area)) +
  geom_histogram(binwidth = 0.0002, fill = "orange", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = median_value), color = "black", linetype = "dashed", size = 1,
             show.legend = TRUE) +
  labs(title = "",
       x = expression("LA " (m^2)), 
       y = "Frequency") +
  annotate("text", x = median_value + 0.002, y = 2300, label = paste("Median:", round(median_value, 6)),
           vjust = -0.5, color = "black", size = 5) +
  # Set x-axis limits from 0 to 500
  xlim(0, 0.01) +
  theme_bw(base_size = 15) 

#-------------------------------
# GT COMPARISON
manual_measurements <- as.data.frame(fread("data/GT_comparison.csv"))
model <- lm(manual_measurements$width_pixels~manual_measurements$pixel_distance)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

gt_comparison_scatter_plot <- ggplot(manual_measurements, aes(x = width_pixels, y = pixel_distance)) +
  geom_point(aes(color = pixel_distance), size = 3, alpha = 0.7) +
  scale_color_viridis_c(option = "C", end = 0.85) +
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  geom_abline(
    intercept = 0, 
    slope = 1, 
    color = "black", linetype = "solid", linewidth = 1
  ) +
  annotate(
    "text", x = min(manual_measurements$width_pixels) , y = max(manual_measurements$pixel_distance) * 0.9, 
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "Petiole width (LM2)",
    y = "Petiole width (manual)",
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------------
merged_dataset$petiole_width <- log(merged_dataset$petiole_width)
merged_dataset$area <- log(merged_dataset$area)
model <- lm(merged_dataset$petiole_width~merged_dataset$area)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

pw_la_comparison_scatter_plot <- ggplot(merged_dataset, aes(x = area, y = petiole_width)) +
  geom_point(aes(color = petiole_width), size = 3, alpha = 0.7) +
  scale_color_viridis_c(option = "C", end = 0.85) +
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(merged_dataset$area) , y = max(merged_dataset$petiole_width) * 0.9, 
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log(LA)",
    y = "log(PW)",
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

#-------------------------------
pdf("plot_for_grant.pdf" ,height=5,width=14)
grid.arrange(pw_la_comparison_scatter_plot, gt_comparison_scatter_plot, ncol=2, nrow = 1)
dev.off()

pdf("plot_for_grant.pdf" ,height=5,width=14)
grid.arrange(petiole_width_hist, LMA_hist, LA_hist, ncol=3, nrow = 1)
dev.off()
