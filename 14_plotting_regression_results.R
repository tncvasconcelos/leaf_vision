library(dplyr)
library(MuMIn)
library(ggplot2)
setwd("~/leaf_vision/")

####### coeficient visualization
# 07a
full_model <- readRDS(file = "models/full_model_07b.rds")
model_set <- readRDS(file = "models/model_set_07b.rds")

## MODEL AVG REGRESSION
# Model averaging of top models within 2 AICc units
model_avg <- model.avg(model_set, subset = delta < 2)

# Extract coefficients and confidence intervals from the model
coef_data <- as.data.frame(summary(model_avg)$coefmat.full)
coef_data <- coef_data[-1,]  # Remove intercept if needed
coef_data$Term <- rownames(coef_data)
# Extract and reorder importance values to match the coefficient terms
importance_values <- sw(model_avg)
importance_values <- importance_values[coef_data$Term]
coef_data$importance_values <- importance_values

coef_data$Term <- c("Aridity Index",
                             "Diurnal Temp. Range",
                             "Max Temp. Warmest Month",
                             "Solar Radiation",
                             "Wind Speed",
                             "Prec. Coldest Month",
                             "Prec. Seasonality",
                             "Prec. Warmest Month",
                             "Elevation",
                             "bio_3")

colnames(coef_data)[4] <-"p_value"
colnames(coef_data)[2] <-"Std_Error"

#
#coef_data <- coef_data[ , c("Estimate", "Std. Error")]
coef_data$Lower <- coef_data$Estimate - 1.96 * coef_data$Std_Error
coef_data$Upper <- coef_data$Estimate + 1.96 * coef_data$Std_Error

# Define a scaling factor to increase the visibility of small estimates
# scaling_factor <- 1e5
# coef_data$Estimate <- coef_data$Estimate * scaling_factor
# coef_data$Lower <- coef_data$Lower * scaling_factor
# coef_data$Upper <- coef_data$Upper * scaling_factor

# Filter out the term with the highest p-value and reorder by p-value in ascending order
coef_data_filtered <- coef_data %>%
  arrange(p_value) %>%
  #filter(p_value < max(p_value)) %>%
  mutate(
    Term = factor(Term, levels = rev(Term)),  # Reorder so lowest p-value is at the top
    color_intensity = ifelse(p_value < 0.05, "high", "low"),  # Adjust color intensity
    color = ifelse(Estimate > 0, "positive", "negative")      # Color by sign of estimate
  )

# Plot
coef_plot <- ggplot(coef_data_filtered, aes(x = Term, y = Estimate, color = color, alpha = color_intensity)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(title = "",
       x = "", y = "Effect Size") +
  scale_color_manual(values = c("positive" = "darkred", "negative" = "darkblue")) +
  scale_alpha_manual(values = c("high" = 1, "low" = 0.5)) +  # Set transparency for non-significant terms
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()  + # Flip coordinates
  guides(color = "none", alpha = "none")  # Remove legend for color and alpha


# Importance values bar plot
importance_plot <- ggplot(coef_data_filtered, aes(x = Term, y = importance_values)) +
  geom_bar(stat = "identity", fill = "grey60") +
  labs(title = "", x = "", y = "Importance Value") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()  # Flip coordinates to match orientation

pdf("FIGURES/plot_for_figure5.pdf" ,height=4,width=8)
grid.arrange(coef_plot, importance_plot, ncol=2, nrow = 1)
dev.off()



