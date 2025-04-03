# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
# setwd("~/leaf_vision/")
library(gridExtra)
library(data.table)
library(ggplot2)

#---------------------------------------
merged_dataset <- read.csv("merged_dataset_final.csv")
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$area))
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$petiole_width))

# gt petiole comparison
manual_measurements <- as.data.frame(fread("data/GT_comparison.csv"))
model <- lm(manual_measurements$width_pixels~manual_measurements$pixel_distance)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

#-------------------------------
# BOOTSTRAP
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
    # Royer equation is log scale
    log_lma <- log10(base_lma)
    # Add error in log space and convert back
    lma_simulations[i] <- 10^(log_lma + rnorm(1, 0, sd_royer))
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
# Add error from Royer equation (R² = 0.55)
variance_explained <- 0.55
sd_royer <- sqrt(var(log(merged_dataset$LMA)) * (1 - variance_explained))

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
}, mc.cores = 4)

# saveRDS(results_list, file = "intermediate_datasets/ci_results_list.RDS")
results_list <- readRDS("ci_results_list.RDS")

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
  geom_point(size = 0.01, color = "darkblue") +
  # Improve labels
  labs(title = "LMA Estimates with 95% Confidence Intervals",
    x = "Leaf (ordered by LMA)", 
    y = "LMA (g/m²)")


#### SOURCES OF VARIATION

# Function to decompose sources of uncertainty
analyze_uncertainty_sources <- function(specimen_data, n_simulations = 1000) {
  # Set up arrays to store results from different scenarios
  lma_full_error <- numeric(n_simulations)      # Both error sources
  lma_width_error_only <- numeric(n_simulations)  # Only width calibration error
  lma_royer_error_only <- numeric(n_simulations)  # Only Royer equation error
  
  for(i in 1:n_simulations) {
    # Base calculation - no errors
    width_base <- coef_model[1] + coef_model[2] * specimen_data$width_pixels
    width_converted_base <- width_base / specimen_data$conversion_mean
    lma_base <- make_LMA(specimen_data$area, width_converted_base)
    log_lma_base <- log(lma_base)
    
    # Width measurement error
    prediction_error_sd <- sqrt(mean(residuals(model)^2))
    width_with_error <- width_base + rnorm(1, 0, prediction_error_sd)
    width_converted_with_error <- width_with_error / specimen_data$conversion_mean
    
    # Scenario 1: Both errors
    lma_full_error[i] <- exp(log(make_LMA(specimen_data$area, width_converted_with_error)) + 
        rnorm(1, 0, sd_royer))
    
    # Scenario 2: Width error only
    lma_width_error_only[i] <- make_LMA(specimen_data$area, width_converted_with_error)
    
    # Scenario 3: Royer error only
    lma_royer_error_only[i] <- exp(log_lma_base + rnorm(1, 0, sd_royer))
  }
  
  # Calculate variances for each scenario
  var_full <- var(lma_full_error)
  var_width <- var(lma_width_error_only)
  var_royer <- var(lma_royer_error_only)
  
  # Calculate confidence interval widths
  ci_width_full <- quantile(lma_full_error, 0.975) - quantile(lma_full_error, 0.025)
  ci_width_width_only <- quantile(lma_width_error_only, 0.975) - quantile(lma_width_error_only, 0.025)
  ci_width_royer_only <- quantile(lma_royer_error_only, 0.975) - quantile(lma_royer_error_only, 0.025)
  
  # Calculate proportion of variance contributed by each source
  prop_width <- var_width / var_full
  prop_royer <- var_royer / var_full
  
  # Return results
  return(list(
    point_estimate = make_LMA(specimen_data$area, width_converted_base),
    variance_total = var_full,
    variance_width = var_width,
    variance_royer = var_royer,
    proportion_width = prop_width,
    proportion_royer = prop_royer,
    ci_width_full = ci_width_full,
    ci_width_width_only = ci_width_width_only,
    ci_width_royer_only = ci_width_royer_only
  ))
}

# Apply to a sample of specimens (or all of them)
sample_size <- min(20, nrow(merged_dataset))
sample_indices <- sample(which(merged_dataset$conversion_mean > 0), sample_size)

uncertainty_results <- lapply(sample_indices, function(i) {
  analyze_uncertainty_sources(merged_dataset[i,])
})

# Summarize results
summary_df <- data.frame(
  specimen_id = sample_indices,
  lma_estimate = sapply(uncertainty_results, function(x) x$point_estimate),
  prop_width_error = sapply(uncertainty_results, function(x) x$proportion_width),
  prop_royer_error = sapply(uncertainty_results, function(x) x$proportion_royer),
  ci_width_full = sapply(uncertainty_results, function(x) x$ci_width_full),
  ci_width_width_only = sapply(uncertainty_results, function(x) x$ci_width_width_only),
  ci_width_royer_only = sapply(uncertainty_results, function(x) x$ci_width_royer_only)
)

# Overall summary
cat("Mean proportion of variance from width measurement:", 
  mean(summary_df$prop_width_error), "\n")
cat("Mean proportion of variance from Royer equation:", 
  mean(summary_df$prop_royer_error), "\n")

# Visualize the contribution of each error source
library(ggplot2)

# Create a stacked bar chart showing relative contributions
contributions_long <- data.frame(
  specimen = rep(summary_df$specimen_id, 2),
  source = c(rep("Width Calibration", nrow(summary_df)), 
    rep("Royer Equation", nrow(summary_df))),
  proportion = c(summary_df$prop_width_error, 
    summary_df$prop_royer_error)
)

p_contributions <- ggplot(contributions_long, aes(x = factor(specimen), y = proportion, fill = source)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Relative Contribution to Uncertainty by Source",
    x = "Specimen ID",
    y = "Proportion of Total Variance") +
  scale_fill_brewer(palette = "Set1")

# Create a boxplot comparing CI widths from different sources
ci_widths_long <- data.frame(
  specimen = rep(summary_df$specimen_id, 3),
  source = c(rep("Both Sources", nrow(summary_df)),
    rep("Width Only", nrow(summary_df)),
    rep("Royer Only", nrow(summary_df))),
  ci_width = c(summary_df$ci_width_full,
    summary_df$ci_width_width_only,
    summary_df$ci_width_royer_only)
)

p_ci_widths <- ggplot(ci_widths_long, aes(x = source, y = ci_width, fill = source)) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "95% CI Width by Error Source",
    x = "",
    y = "Confidence Interval Width (g/m²)") +
  scale_fill_brewer(palette = "Set1")

# Display plots
print(p_contributions)
print(p_ci_widths)
