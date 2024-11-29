#------------------------------------
# Load necessary libraries
library(phylolm)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)

load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))


#--------------------------
# BIO2
model <- phylolm(lma~bio_2, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_2_plot <- ggplot(dat, aes(x = bio_2, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_2), y = 6.9,
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Mean Diurnal Range",
    y = expression("log LMA " (g/m^2))
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

#--------------------------
# BIO3
model <- phylolm(lma~bio_3, data=dat, phy=phy, model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_3_plot <- ggplot(dat, aes(x = bio_3, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_2) , y = 6.9,
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Isothermality",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
#-------------------------
# BIO5
model <- phylolm(lma~bio_5, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_5_plot <- ggplot(dat, aes(x = bio_5, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_5) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Max Temperature of the Warmest Month ",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
# BIO15
model <- phylolm(lma~bio_15, data=dat, phy=phy, model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_15_plot <- ggplot(dat, aes(x = bio_15, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_5) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Precipitation Seasonality ",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
# BIO18
model <- phylolm(lma~bio_18, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_18_plot <- ggplot(dat, aes(x = bio_18, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_5) , y =6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Precipitation of Warmest Quarter",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
# BIO19
model <- phylolm(lma~bio_19, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_19_plot <- ggplot(dat, aes(x = bio_19, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_5) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Precipitation of Coldest Quarter",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
#-------------------------
# AI
model <- phylolm(lma~ai, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_ai_plot <- ggplot(dat, aes(x = ai, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$ai) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Aridity Index",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
#-------------------------
# S-RAD
model <- phylolm(lma~srad, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_srad_plot <- ggplot(dat, aes(x = srad, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$bio_5) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Solar Radiation",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
#-------------------------
# Wind
model <- phylolm(lma~wind, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_wind_plot <- ggplot(dat, aes(x = wind, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$wind) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Mean Annual Wind speed",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
#-------------------------
# Alt
model <- phylolm(lma~alt, data=dat, phy=phy,model = "lambda", REML = FALSE)
coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_alt_plot <- ggplot(dat, aes(x = alt, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1], 
    slope = coef(model)[2], 
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  # annotate(
  #   "text", x = min(dat$wind) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  labs(
    x = "Elevation (m)",
    y = expression("log LMA " (g/m^2)),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")
#-------------------------
pdf("FIGURES/plot_for_figure4a.pdf" ,height=12,width=15)
grid.arrange(lma_bio_2_plot,
             lma_bio_5_plot,
             lma_bio_15_plot,
             lma_bio_18_plot,
             lma_bio_19_plot,
             lma_srad_plot,
             lma_ai_plot, 
             lma_wind_plot,
             lma_alt_plot,
             ncol=3, nrow = 3)
dev.off()
###

#lma_bio_3_plot

