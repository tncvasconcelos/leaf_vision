#------------------------------------
# Load necessary libraries
library(phylolm)
library(phytools)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)
library(lmtest)
# rm(list=ls())

merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset$LMA <- log((exp(merged_dataset$LMA) * 100))

merged_dataset <- subset(merged_dataset, merged_dataset$lat<40 & merged_dataset$lat>-40)


merged_dataset$lat <- abs(merged_dataset$lat)

# plot(merged_dataset$LMA~log(abs(merged_dataset$lat)))
model <- lm(merged_dataset$LMA~merged_dataset$lat)
# abline(model, col="red")

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_lat_plot <- ggplot(merged_dataset, aes(x = lat, y = LMA)) +
  geom_point(aes(color = LMA), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "C", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(merged_dataset$lat) , y = max(merged_dataset$LMA),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log Latitude",
    y = "log LMA",
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# setwd("~/leaf_vision/")
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
phy <- multi2di(phy)
phy <- force.ultrametric(phy)

#--------------------------
#--------------------------
# BIO1
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
model <- phylolm(lma~bio_1, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_1_plot <- ggplot(dat, aes(x = bio_1, y = lma)) +
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = bio_1 - bio_1_se, xmax = bio_1 + bio_1_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_1) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "MAT",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$bio_1_se <- log(dat$bio_1_se) 
dat <- subset(dat, dat$bio_1_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~bio_1_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_1_se_plot <- ggplot(dat, aes(x = bio_1_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_1_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log MAT SE",
    y = expression("log LMA SE"),
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
# BIO12
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
model <- phylolm(lma~bio_12, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_12_plot <- ggplot(dat, aes(x = bio_12, y = lma)) +
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = bio_12 - bio_12_se, xmax = bio_12 + bio_12_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_12) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "MAP",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$bio_12_se <- log(dat$bio_12_se) 
dat <- subset(dat, dat$bio_12_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~bio_12_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_12_se_plot <- ggplot(dat, aes(x = bio_12_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_12_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log MAP SE",
    y = expression("log LMA SE"),
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
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
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
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = bio_15 - bio_15_se, xmax = bio_15 + bio_15_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_15) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "Prec. Seasonality",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$bio_15_se <- log(dat$bio_15_se) 
dat <- subset(dat, dat$bio_15_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~bio_15_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_15_se_plot <- ggplot(dat, aes(x = bio_15_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_15_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log Prec. Seasonality SE",
    y = expression("log LMA SE"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

#-------------------------
# BIO4
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
model <- phylolm(lma~bio_4, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_4_plot <- ggplot(dat, aes(x = bio_4, y = lma)) +
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = bio_4 - bio_4_se, xmax = bio_4 + bio_4_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_4) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "Temp. Seasonality",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$bio_4_se <- log(dat$bio_4_se) 
dat <- subset(dat, dat$bio_4_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~bio_4_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_4_se_plot <- ggplot(dat, aes(x = bio_4_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$bio_4_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log Temp. Seasonality SE",
    y = expression("log LMA SE"),
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
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
model <- phylolm(lma~ai, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_ai_plot <- ggplot(dat, aes(x = ai, y = lma)) +
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = ai - ai_se, xmax = ai + ai_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$ai) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "Aridity Index",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$ai_se <- log(dat$ai_se) 
dat <- subset(dat, dat$ai_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~ai_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_ai_se_plot <- ggplot(dat, aes(x = ai_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$ai_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log Aridity Index SE",
    y = expression("log LMA SE"),
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
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
model <- phylolm(lma~srad, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_srad_plot <- ggplot(dat, aes(x = srad, y = lma)) +
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = srad - srad_se, xmax = srad + srad_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$srad) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "Solar Radiation",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$srad_se <- log(dat$srad_se) 
dat <- subset(dat, dat$srad_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~srad_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_srad_se_plot <- ggplot(dat, aes(x = srad_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$srad_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log Solar Radiation SE",
    y = expression("log LMA SE"),
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
load("results/data_subset_for_plots.Rsave")
dat$lma <- log((exp(dat$lma) * 100))
# species means plot
model <- phylolm(lma~wind, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_wind_plot <- ggplot(dat, aes(x = wind, y = lma)) +
  geom_errorbar(aes(ymin = lma - se, ymax = lma + se),
                width = 0, color = "black", alpha = 0.3) +   # Vertical error bars
  geom_errorbarh(aes(xmin = wind - wind_se, xmax = wind + wind_se),
                 height = 0, color = "black", alpha = 0.3) + # Horizontal error bars
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$wind) , y = 6.9,
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "Wind Speed",
    y = expression("log LMA"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

# standart error plot
dat$wind_se <- log(dat$wind_se) 
dat <- subset(dat, dat$wind_se!=-Inf)
dat$se <- log(dat$se) 

model <- phylolm(se~wind_se, data=dat, phy=phy, model = "lambda", REML = FALSE)

coef_model <- coef(model)  # Intercept and slope
r2 <- summary(model)$r.squared  # R-squared
p_value <- summary(model)$coefficients[2, 4]  # p-value for the slope

# Create a text label with R² and p-value
label_text <- paste0(
  "R² = ", round(r2, 3), 
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_wind_se_plot <- ggplot(dat, aes(x = wind_se, y = se)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  scale_color_viridis_c(option = "B", end = 0.85) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  # Add a custom linear model trend line using the equation from `model`
  geom_abline(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text", x = min(dat$wind_se) , y = max(dat$se),
    label = label_text, size = 5, hjust = 0, vjust = 1
  ) +
  labs(
    x = "log Wind Speed SE",
    y = expression("log LMA SE"),
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")


#-------------------------
pdf("FIGURES/plot_for_figure_SI_climatic_vars1.pdf" ,height=10,width=10)
grid.arrange(lma_bio_1_plot, lma_bio_1_se_plot,
             lma_bio_4_plot, lma_bio_4_se_plot,
             ncol=2, nrow = 2,
             widths = c(1,1))
dev.off()

pdf("FIGURES/plot_for_figure_SI_climatic_vars2.pdf" ,height=10,width=10)
grid.arrange(lma_bio_12_plot, lma_bio_12_se_plot,
             lma_bio_15_plot, lma_bio_15_se_plot,
             ncol=2, nrow = 2,
             widths = c(1,1))
dev.off()

pdf("FIGURES/plot_for_figure_SI_climatic_vars3.pdf" ,height=15,width=10)
grid.arrange(lma_ai_plot, lma_ai_se_plot,
             lma_srad_plot, lma_srad_se_plot,
             lma_wind_plot, lma_wind_se_plot,
             ncol=2, nrow =3,
             widths = c(1,1))
dev.off()


