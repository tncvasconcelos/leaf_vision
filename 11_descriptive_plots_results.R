# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
# setwd("~/leaf_vision/")
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
# model <- lm(log(merged_dataset$area)~log(merged_dataset$petiole_width))
# summary(model)
# plot(log(merged_dataset$area)~log(merged_dataset$petiole_width), 
#      xlab="log(petiole width cm)", ylab="log(leaf area cm^2)")
# 
# abline(model, col="red")
#dev.off()


#-------------------------------
# FIGURE 3
#-------------------------------
# gt petiole comparison
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


#-------------------------------
# FIGURE BOOTSTRAP
#-------------------------------
# Function to create confidence intervals for individual LMA estimates
calculate_individual_lma_with_ci <- function(specimen_data, n_simulations = 1000) {
  # Initialize storage for simulated LMA values
  lma_simulations <- numeric(n_simulations)
  
  for(i in 1:n_simulations) {
    # Step 1: Apply calibration model with error
    calibrated_width <- coef_model[1] + coef_model[2] * specimen_data$width_pixels
    prediction_error_sd <- sqrt(mean(residuals(model)^2))
    calibrated_width <- calibrated_width + rnorm(1, 0, prediction_error_sd)
    
    # Convert from pixels to meters
    converted_width <- calibrated_width / specimen_data$conversion_mean
    
    # Step 2: Apply Royer equation
    base_lma <- make_LMA(specimen_data$area, converted_width)
    log_lma <- log(base_lma)
    
    # Add error from Royer equation (R² = 0.55)
    variance_explained <- 0.55
    # Standard error of the Royer regression (can be approximated)
    # Using the coefficient from Royer equation
    se_royer <- 0.382 * sqrt((1 - variance_explained) / variance_explained)  
    
    # Add error in log space and convert back
    lma_simulations[i] <- exp(log_lma + rnorm(1, 0, se_royer))
  }
  
  # Calculate point estimate and confidence interval
  point_estimate <- make_LMA(
    specimen_data$area, 
    (coef_model[1] + coef_model[2] * specimen_data$width_pixels) / specimen_data$conversion_mean
  )
  
  ci_lower <- quantile(lma_simulations, 0.025)
  ci_upper <- quantile(lma_simulations, 0.975)
  
  return(list(
    estimate = point_estimate,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    simulations = lma_simulations
  ))
}

make_LMA <- function(leaf_area, petiole_width) {
  LMA = 3.07 + 0.382 * log(petiole_width^2 / leaf_area)
  return(exp(LMA))
}

library(parallel)
# Apply to all specimens in your dataset
results_list <- mclapply(1:nrow(merged_dataset), function(i) {
  # Skip rows with conversion_mean of 0
  if(merged_dataset$conversion_mean[i] == 0) {
    return(list(
      estimate = NA,
      lower_ci = NA,
      upper_ci = NA,
      simulations = NA
    ))
  }
  
  calculate_individual_lma_with_ci(merged_dataset[i,])
}, mc.cores = 10)

# Extract results into a data frame
output <- data.frame(
  specimen_id = merged_dataset$component_name, 
  lma_estimate = sapply(results_list, function(x) x$estimate),
  lma_lower_ci = sapply(results_list, function(x) x$lower_ci),
  lma_upper_ci = sapply(results_list, function(x) x$upper_ci)
)

# Filter out any rows with NAs (from rows with conversion_mean = 0)
output <- output[!is.na(output$lma_estimate),]

# Print some summary statistics
cat("Mean LMA across all specimens:", mean(output$lma_estimate), "g/m²\n")
cat("Mean CI width:", mean(output$lma_upper_ci - output$lma_lower_ci), "g/m²\n")

# Create histogram of LMA estimates with error bars
library(ggplot2)
# Sort by LMA value for better visualization
output_sorted <- output[order(output$lma_estimate, decreasing = TRUE),]
output_sorted$order <- 1:nrow(output_sorted)

# Create the plot with semi-transparent error bars
ggplot(output_sorted, aes(x = order, y = lma_estimate)) +
  # Add a background grid only for the y-axis for better readability
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  # Add the error bars first (so points appear on top)
  geom_errorbar(aes(ymin = lma_lower_ci, ymax = lma_upper_ci), 
    width = 0.2, 
    color = "steelblue",  # Less intense blue
    alpha = 0.5) +
  # Add points last so they're on top
  geom_point(size = 1, color = "darkblue") +
  # Improve labels
  labs(title = "LMA Estimates with 95% Confidence Intervals",
    x = "Specimens (ordered by LMA)", 
    y = "LMA (g/m²)")

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
# load("results/data_subset_w_leaf_phenology.Rsave")
# data_subset$lma <- exp(data_subset$lma)*100
# 
# median(data_subset$lma[data_subset$deciduousness=="evergreen"])
# median(data_subset$lma[data_subset$deciduousness=="deciduous"])
# 
# 
# pdf("FIGURES/boxplot_deciduousness.pdf", height=3, width=6)
# ggplot(data_subset, aes(x = deciduousness, y = lma, fill = deciduousness)) +
#   geom_boxplot(outlier.shape = NA, width = 0.6) +  # Box plot without outliers
#   geom_jitter(aes(color = deciduousness), width = 0.2, size = 1.5, alpha = 0.3) +  # Jitter points with matching color and increased transparency
#   theme_bw() +
#   labs(x = "Deciduousness", y = "LMA") +
#   scale_fill_manual(values = c("evergreen" = "darkgreen", "deciduous" = "orange")) +  # Custom fill colors
#   scale_color_manual(values = c("evergreen" = "darkgreen", "deciduous" = "orange")) +  # Matching color for jitter points
#   theme(legend.position = "none")  # Remove legend if not needed
# dev.off()
