library(ape)
library(MuMIn)
library(rr2)
library(phylolm)
setwd("~/leaf_vision/")

####### coeficient visualization
# 07a
full_model <- readRDS(file = "models/full_model_07a.rds")
model_set <- readRDS(file = "models/model_set_07a.rds")
phy <- read.tree("trees/phy_070a.tre")
data_subset <- dat <- full_model$model.frame
full_formula <- full_model$formula

# Extract coefficients and confidence intervals from the model
coef_data <- as.data.frame(summary(model_avg)$coefmat.full)
coef_data <- coef_data[ , c("Estimate", "Std. Error")]
coef_data$Lower <- coef_data$Estimate - 1.96 * coef_data$`Std. Error`
coef_data$Upper <- coef_data$Estimate + 1.96 * coef_data$`Std. Error`
coef_data$Term <- rownames(coef_data)
coef_data <- coef_data[-1,]  # Remove intercept if needed

# Extract and reorder importance values to match the coefficient terms
importance_values <- sw(model_avg)
importance_values <- importance_values[coef_data$Term]

# Set up layout for side-by-side plots
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

# Left Panel: Model-Averaged Coefficients with Confidence Intervals
plot(coef_data$Estimate, 1:nrow(coef_data), xlim = range(coef_data$Lower, coef_data$Upper),
  pch = 19, xlab = "Estimate", ylab = "", yaxt = "n", main = "Model-Averaged Coefficients")
axis(2, at = 1:nrow(coef_data), labels = coef_data$Term, las = 1)
segments(coef_data$Lower, 1:nrow(coef_data), coef_data$Upper, 1:nrow(coef_data))
abline(v = 0, lty = 2, col = "grey")

# Right Panel: Predictor Importance without Labels, matching the y-axis order
barplot(importance_values, horiz = TRUE, las = 1,
  xlab = "Per-variable\nsum of model weights", col = "forestgreen", main = "Predictor Importance",
  names.arg = NULL, space = 0.9, )  # Adjust spacing to match alignment


#### many regression lines
# Extract model-averaged coefficients
coef_data <- as.data.frame(summary(model_avg)$coefmat.full)
coef_data <- coef_data[ , "Estimate"]
names(coef_data) <- rownames(summary(model_avg)$coefmat.full)

# Assuming your original dataset is `data_subset` and dependent variable is `lma`
dependent_var <- "lma"
predictors <- names(coef_data)[-1]  # Exclude intercept

# Set up multi-panel layout
par(mfrow = c(2, ceiling(length(predictors) / 2)), mar = c(4, 4, 2, 1))

for (pred in predictors) {
  # Plot raw values of predictor vs. dependent variable
  plot(data_subset[[pred]], data_subset[[dependent_var]], main = paste("Effect of", pred),
    xlab = pred, ylab = dependent_var, pch = 19, col = "blue")
  
  # Add the model-averaged line using the model-averaged intercept and slope for each predictor
  intercept <- coef_data["(Intercept)"]
  slope <- coef_data[pred]
  abline(a = intercept, b = slope, col = "red", lwd = 2)
}
