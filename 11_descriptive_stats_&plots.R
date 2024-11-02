# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
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
    x = "LM2 Petiole width (pixels)",
    y = "manual Petiole width (pixels)",
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# #-------------------------------
# merged_dataset$petiole_width <- log(merged_dataset$petiole_width)
# merged_dataset$area <- log(merged_dataset$area)
# model <- lm(merged_dataset$petiole_width~merged_dataset$area)
# coef_model <- coef(model)  # Intercept and slope
# r2 <- summary(model)$r.squared  # R-squared
# p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope
# 
# # Create a text label with R² and p-value
# label_text <- paste0(
#   "R² = ", round(r2, 3), 
#   "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
# )
# 
# pw_la_comparison_scatter_plot <- ggplot(merged_dataset, aes(x = area, y = petiole_width)) +
#   geom_point(aes(color = petiole_width), size = 3, alpha = 0.7) +
#   scale_color_viridis_c(option = "C", end = 0.85) +
#   # Add a custom linear model trend line using the equation from `model`
#   geom_abline(
#     intercept = coef(model)[1], 
#     slope = coef(model)[2], 
#     color = "black", linetype = "dashed", linewidth = 1
#   ) +
#   annotate(
#     "text", x = min(merged_dataset$area) , y = max(merged_dataset$petiole_width) * 0.9, 
#     label = label_text, size = 5, hjust = 0, vjust = 1
#   ) +
#   labs(
#     x = "log(LA)",
#     y = "log(PW)",
#   ) +
#   theme_bw(base_size = 15) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     panel.grid.major = element_line(color = "grey85"),
#     legend.position = ""
#   ) +
#   ggtitle("")

#-------------------------------
pdf("FIGURES/PW_comparison.pdf" ,height=5,width=5)
gt_comparison_scatter_plot
dev.off()






# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
library(gridExtra)
library(data.table)
library(ggplot2)

#---------------------------------------
merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$area))
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$petiole_width))

#-------------------------------
# Histogram for PW
median_value1 <- median(merged_dataset$petiole_width)
petiole_width_hist <- ggplot(merged_dataset, aes(x = petiole_width)) +
  geom_histogram(binwidth = 0.018, fill = "orange", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = median_value1), color = "black", linetype = "dashed", linewidth = 1,
             show.legend = TRUE) +
  labs(title = "",
       x = "LM2 Petiole width (cm)", 
       y = "Frequency") +
  annotate("text", x = median_value1 + 0.15, y = 3000, label = paste("Median:", round(median_value1, 2)),
           vjust = -0.5, color = "black", size = 5) +
  # Set x-axis limits from 0 to 500
  theme_bw(base_size = 15) 

#-------------------------------
# Histogram for LMA
merged_dataset$LMA_for_plot <- merged_dataset$LMA*100
median_value2 <- median(merged_dataset$LMA_for_plot)

LMA_hist <- ggplot(merged_dataset, aes(x = LMA_for_plot)) +
  geom_histogram(binwidth = 10, fill = "darkorchid", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = median_value2), color = "black", linetype = "dashed", size = 1,
             show.legend = TRUE) +
  labs(title = "",
       x = expression("LMA " (g/m^2)), 
       y = "Frequency") +
  annotate("text", x = median_value2 + 100, y = 1570, label = paste("Median:", round(median_value2, 1)),
           vjust = -0.5, color = "black", size = 5) +
  # Set x-axis limits from 0 to 500
  xlim(0, 500) +
  theme_bw(base_size = 15) 
#-------------------------------
# Histogram for area
median_value3 <- median(merged_dataset$area)
LA_hist <- ggplot(merged_dataset, aes(x = area)) +
  geom_histogram(binwidth = 3, fill = "red3", color = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = median_value3), color = "black", linetype = "dashed", size = 1,
             show.legend = TRUE) +
  labs(title = "",
       x = expression("LM2 Leaf Area " (cm^2)), 
       y = "Frequency") +
  annotate("text", x = median_value3 + 35, y = 4850, label = paste("Median:", round(median_value3, 2)),
           vjust = -0.5, color = "black", size = 5) +
  # Set x-axis limits from 0 to 500
  xlim(0, 150) +
  theme_bw(base_size = 15) 

pdf("FIGURES/plot_for_figure3_test.pdf" ,height=10,width=10)
grid.arrange(petiole_width_hist, LA_hist, gt_comparison_scatter_plot, LMA_hist, ncol=2, nrow = 2)
dev.off()

#-------------------------------
