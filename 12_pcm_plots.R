


library(ggplot2)


# Split data by binary variable
df0 <- subset(df, binary_var == 0)
df1 <- subset(df, binary_var == 1)

# Calculate regression lines (replace with your custom model)
# Example with hypothetical custom_model function:
model0 <- my_custom_model(y ~ x, data = df0)  # Fit model for binary_var = 0
model1 <- my_custom_model(y ~ x, data = df1)  # Fit model for binary_var = 1

# Predict values based on models
df0$y_pred <- predict(model0, newdata = df0)
df1$y_pred <- predict(model1, newdata = df1)

# Combine predictions for plotting
df_pred <- rbind(df0, df1)

# Plot with ggplot2
ggplot(df, aes(x = x, y = y, color = factor(binary_var))) +
  geom_point() +  # Scatter plot points
  geom_line(data = df_pred, aes(x = x, y = y_pred, color = factor(binary_var))) +  # Custom regression lines
  labs(color = "Binary Variable") +
  theme_minimal()

ggplot(data_subset, aes(x = srad, y = lma, color = factor(deciduousness))) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", se = FALSE) +  # Regression lines without confidence intervals
  labs(color = "Binary Variable") +  # Label for the legend
  theme_minimal()  # Minimal theme for better aesthetics

