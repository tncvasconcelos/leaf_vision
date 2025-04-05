#------------------------------------
# Load necessary libraries
library(phylolm)
library(phytools)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)
library(lmtest)
library(geomorph)

# setwd("~/leaf_vision/")
# dat <- read.csv("data/merged_dataset_final.csv")
# load("results/data_subset_for_plots.Rsave")
# dat$lma <- log((exp(dat$lma) * 100))
# phy <- multi2di(phy)
# phy <- force.ultrametric(phy)

merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset$genus_species <- gsub(" ", "_", merged_dataset$genus_species)
merged_dataset <- merged_dataset[!is.na(merged_dataset$genus_species),]
tre <- read.tree("trees/GBMB.tre")
missing_sp <- unique(merged_dataset$genus_species)[!unique(merged_dataset$genus_species) %in% tre$tip.label]
dat <- merged_dataset[!merged_dataset$genus_species %in% missing_sp,]
phy <- keep.tip(tre, tre$tip.label[tre$tip.label %in% dat$genus_species])
phy <- force.ultrametric(phy)
phy$node.label <- NULL

#--------------------------
# BIO1
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$bio_1), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
bio1_dat <- data.frame(sp = tmp$Group.2, bio_1 = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ bio_1,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = bio1_dat,
  print.progress = TRUE
)
save(model, file="bio1_extd_model.Rsave")

load(file="bio1_extd_model.Rsave")
lrtest_bio1 <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_bio1$table$`Pr(>F)`[1]

# Create a text label with R² and p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_1_plot <- ggplot(bio1_dat, aes(x = bio_1, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Mean Annual Temperature",
    y = expression("log LMA " (g/m^2))
  ) +
  # annotate(
  #   "text", x = min(dat$wind) , y = 6.9, 
  #   label = label_text, size = 5, hjust = 0, vjust = 1
  # ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")


#--------------------------
# BIO12
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$bio_12), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
bio12_dat <- data.frame(sp = tmp$Group.2, bio_12 = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ bio_12,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = bio12_dat,
  print.progress = TRUE
)
save(model, file="bio12_extd_model.Rsave")

load(file="bio12_extd_model.Rsave")
lrtest_bio12 <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_bio12$table$`Pr(>F)`[1]

# Create a text label with p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_12_plot <- ggplot(bio12_dat, aes(x = bio_12, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Annual Precipitation",
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
# Solar radiation
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$srad), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
srad_dat <- data.frame(sp = tmp$Group.2, srad = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ srad,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = srad_dat,
  print.progress = TRUE
)
save(model, file="srad_extd_model.Rsave")

load(file="srad_extd_model.Rsave")
lrtest_srad <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_srad$table$`Pr(>F)`[1]

# Create a text label with p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_srad_plot <- ggplot(srad_dat, aes(x = srad, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Solar Radiation",
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
# BIO4
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$bio_4), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
bio4_dat <- data.frame(sp = tmp$Group.2, bio_4 = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ bio_4,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = bio4_dat,
  print.progress = TRUE
)
save(model, file="bio4_extd_model.Rsave")

load(file="bio4_extd_model.Rsave")
lrtest_bio4 <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_bio4$table$`Pr(>F)`[1]

# Create a text label with p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_4_plot <- ggplot(bio4_dat, aes(x = bio_4, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Temperature Seasonality",
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
# BIO15
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$bio_15), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
bio15_dat <- data.frame(sp = tmp$Group.2, bio_15 = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ bio_15,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = bio15_dat,
  print.progress = TRUE
)
save(model, file="bio15_extd_model.Rsave")

load(file="bio15_extd_model.Rsave")
lrtest_bio15 <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_bio15$table$`Pr(>F)`[1]

# Create a text label with p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_bio_15_plot <- ggplot(bio15_dat, aes(x = bio_15, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Precipitation Seasonality",
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
# Wind speed
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$wind), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
wind_dat <- data.frame(sp = tmp$Group.2, wind = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ wind,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = wind_dat,
  print.progress = TRUE
)
save(model, file="wind_extd_model.Rsave")

load(file="wind_extd_model.Rsave")
lrtest_wind <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_wind$table$`Pr(>F)`[1]

# Create a text label with p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_wind_plot <- ggplot(wind_dat, aes(x = wind, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Wind Speed",
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
# AI (Aridity Index)
tmp <- data.frame(aggregate(cbind(dat$LMA, dat$ai), list(dat$filename, dat$genus_species),
                            function(x) c(mean(x), mean(log(x)))))
ai_dat <- data.frame(sp = tmp$Group.2, ai = tmp$V2[,1], lma = log(exp(tmp$V1[,2])*100))
model <- extended.pgls(
  lma ~ ai,
  phy = phy,
  Cov = NULL,
  species = "sp",
  delta = 0.001,
  gamma = c("sample", "equal"),
  iter = 999,
  seed = NULL,
  int.first = FALSE,
  turbo = FALSE,
  Parallel = FALSE,
  verbose = FALSE,
  data = ai_dat,
  print.progress = TRUE
)
save(model, file="ai_extd_model.Rsave")

load(file="ai_extd_model.Rsave")
lrtest_ai <- anova(model)
coef_model <- coef(model)
p_value <- lrtest_ai$table$`Pr(>F)`[1]

# Create a text label with p-value
label_text <- paste0(
  "\n", "p = ", format.pval(p_value, digits = 3, eps = 0.001)
)

lma_ai_plot <- ggplot(ai_dat, aes(x = ai, y = lma)) +
  geom_point(aes(color = lma), size = 2, alpha = 0.5) +
  geom_density_2d(color = "black", linewidth = 0.4) +  # Add contour lines
  scale_color_viridis_c(option = "D", end = 0.85) +
  geom_abline(
    intercept = coef_model[1],
    slope = coef_model[2],
    color = "black", linetype = "dashed", linewidth = 1
  ) +
  labs(
    x = "Aridity Index",
    y = expression("log LMA " (g/m^2))
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = ""
  ) +
  ggtitle("")

#------------------------------------
# Now create the final plot with grid.arrange

pdf("plot_for_figure4a.pdf", height=12, width=15)
grid.arrange(lma_bio_1_plot,
             lma_bio_12_plot,
             lma_srad_plot,
             lma_bio_4_plot,
             lma_bio_15_plot,
             lma_wind_plot,
             lma_ai_plot,
             ncol=3, nrow = 3)
dev.off()
