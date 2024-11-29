# standardized values of LMA
library(ape)
library(phylolm)
library(geiger)
library(phytools)
library(scales)
library(RColorBrewer)

setwd("~/leaf_vision/")

merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset$genus_species <- gsub(" ", "_", merged_dataset$genus_species)
tre <- read.tree("trees/GBMB.tre")

lma_results <- aggregate(merged_dataset$LMA, list(merged_dataset$genus_species), 
  FUN = function(x) c(mean(log(x)), sd(log(x))/length(x)))
lma_results <- data.frame(sp = lma_results$Group.1,
  lma = lma_results$x[,1],
  se =  lma_results$x[,2])

merged_dataset <- merged_dataset[!duplicated(merged_dataset$filename),]

merged_dataset_2 <- merged_dataset[,!colnames(merged_dataset) == "super_biome"]
merged_dataset_2 <- merged_dataset_2[,!colnames(merged_dataset_2) == "deciduousness"]
merged_dataset_2 <- merged_dataset_2[,!colnames(merged_dataset_2) == "leaf_phenology"]
focal_cols <- grep("bio_1$", colnames(merged_dataset_2)):ncol(merged_dataset_2)
climate_data <- aggregate(merged_dataset_2[,focal_cols], list(merged_dataset_2$genus_species), 
  FUN = function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm = TRUE)/length(na.omit(x))))
climate_data <- do.call(cbind, lapply(climate_data, function(x) as.data.frame(x)))[,-1]
colnames(climate_data) <- gsub("\\.V1", "", colnames(climate_data))
colnames(climate_data) <- gsub("\\.V2", "_se", colnames(climate_data))
dat <- cbind(lma_results, climate_data)

missing_sp <- dat$sp[!dat$sp %in% tre$tip.label]
dat <- dat[dat$sp %in% tre$tip.label, ]

phy <- keep.tip(tre, dat$sp)
phy <- ladderize(phy)
phy$node.label <- NULL
rownames(dat) <- dat$sp
dat <- dat[phy$tip.label,]
H <- max(node.depth.edgelength(phy))

# List of columns to standardize (excluding species and standard error columns)
vars_to_scale <- grep("_se$", names(dat), invert = TRUE, value = TRUE)
vars_to_scale <- vars_to_scale[!vars_to_scale %in% c("sp", "se")]

# Apply scaling: (value - mean) / standard deviation
dat[vars_to_scale] <- scale(dat[vars_to_scale])

# Check the first few rows to confirm the scaling
head(dat)


write.csv(dat, file = "data/lma_dat_b.csv", row.names = FALSE)

rescale_values <- function(x, a, b) {
  (a + (x - min(x)) / (max(x) - min(x)) * (b - a))
}
plot_tip_values <- setNames(rescale_values(dat$lma, H*1.01, H*1.1), dat$sp)

plot.phylo(phy, show.tip.label = FALSE, no.margin = TRUE, direction = "upwards", 
  y.lim=c(0, H*1.1))
for(i in seq_along(plot_tip_values)){
  segments(i, H*1.01, i, plot_tip_values[i], lwd = 0.5)
}

# Load libraries
library(ape)       
library(nlme)      
library(caper)     
library(phytools)  
library(ggplot2)   
library(corrplot)  

# View the first few rows of the data
head(dat)

# Get a summary of the data
summary(dat)

# Check for missing values
colSums(is.na(dat))

# Histogram of LMA
ggplot(dat, aes(x = lma)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Leaf Mass Area (LMA)", x = "LMA", y = "Frequency")
dev.off()
# Melt the data
library(reshape2)
predictor_vars <- grep("_se$", names(dat), invert = TRUE, value = TRUE)
predictor_vars <- predictor_vars[!predictor_vars %in% c("sp", "lma", "se")]
predictor_vars <- predictor_vars[predictor_vars!="et0"]

# Pairwise scatter plots
# pdf("plots/pair_scatter_plots_b.pdf", height = 14, width = 14)
# pairs(dat[, c("lma", predictor_vars)], main = "Pairwise Scatter Plots")
# dev.off()

## CORRELATION EXAMINATION
# Compute correlation matrix
cor_mat <- cor(dat[, predictor_vars], use = "complete.obs")

# Visualize the correlation matrix
# pdf("plots/corr_matrix_b.pdf")
# corrplot(cor_mat, method = "color", type = "upper", tl.cex = 0.8, tl.col = "black")
# dev.off()

library(caret)
library(phylolm)
# # Find variables with correlation higher than 0.7
# high_cor <- findCorrelation(cor_mat, cutoff = 0.7)
# reduced_vars <- predictor_vars[-high_cor]
# 
# # Updated list of variables
# print(reduced_vars)
# 
# ## PCA EXAMINATION
# # Perform PCA on standardized variables
# pca_res <- prcomp(dat[, predictor_vars], scale. = TRUE)
# 
# # Scree plot to determine the number of principal components to retain
# plot(pca_res, type = "l")
# 
# # Get the proportion of variance explained
# cumsum(summary(pca_res)$importance[2, ])


reduced_vars <- c( "bio_1" , "bio_12" ,"bio_4","bio_15", "wind", "srad" , "ai")

## FULL PHYLO REGRESSION
data_subset <- dat[, c("lma", reduced_vars)]
formula_full <- as.formula(paste("lma ~", paste(reduced_vars, collapse = " + ")))
full_model <- phylolm(formula_full, phy = phy, data = data_subset, model = "lambda", REML = FALSE)

# Summary of the model
summary(full_model)
coef_summary <- summary(full_model)$coefficients
saveRDS(full_model, file = "models/full_model_07b.rds")
write.csv(coef_summary, file = "tables/model_coefficients_07b.csv", row.names = TRUE)


## DATA DREDGE REGRESSION
# Load MuMIn package for model selection
library(MuMIn)

# Perform model selection
model_set <- dredge(full_model, trace = TRUE, rank = "AICc")
write.csv(round(as.data.frame(model_set),3), file="model_set_standartized.csv")
saveRDS(model_set, file = "models/model_set_07b.rds")
model_set <- readRDS(file = "models/model_set_07b.rds")

# View the top models
head(model_set)

## MODEL AVG REGRESSION
# Model averaging of top models within 2 AICc units
avg_model <- model.avg(model_set, subset = delta < 2)

# Summary of the averaged model
summary(avg_model)
model_summary_output <- capture.output(summary(avg_model))
writeLines(model_summary_output, "tables/model_summary_output_07b.csv")

coef_summary <- summary(avg_model)$coefficients
write.csv(avg_model$msTable, file = "tables/model_fits_07b.csv", row.names = TRUE)
write.csv(coef_summary, file = "tables/modelavg_coefficients_07b.csv", row.names = TRUE)

## VALIDATION
residuals_phylolm <- residuals(full_model)
fitted_values_phylolm <- fitted(full_model)
plot(fitted_values_phylolm, residuals_phylolm,
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

### REPEAT ANALYSIS WITH LEAF PHENOLOGY
deciduousness <- data.frame(sp = merged_dataset$genus_species, 
  deciduousness = as.factor(merged_dataset$leaf_phenology))
deciduousness <- deciduousness[!is.na(deciduousness$deciduousness),]
deciduousness_df <- aggregate(deciduousness$deciduousness, by = list(deciduousness$sp), 
  FUN = function(x) (table(factor(x, levels = levels(deciduousness$deciduousness)))))
deciduousness_df <- as.data.frame(do.call(cbind, deciduousness_df))
colnames(deciduousness_df) <- c("sp", levels(deciduousness$deciduousness))
deciduousness_vec <- setNames(as.factor(levels(deciduousness$deciduousness)[apply(deciduousness_df[,-1], 1, which.max)]), deciduousness_df$sp)
deciduousness_vec <- deciduousness_vec[phy$tip.label]
data_subset <- dat[, c("lma", reduced_vars)]
data_subset$deciduousness <- deciduousness_vec
phy <- drop.tip(phy, rownames(data_subset)[which(is.na(data_subset$deciduousness))])
data_subset <- data_subset[!is.na(data_subset$deciduousness), ]
formula_full <- as.formula(paste("lma ~", paste(c(reduced_vars, "deciduousness"), collapse = " + ")))
full_model <- phylolm(formula_full, phy = phy, data = data_subset, model = "lambda", REML = FALSE)

aggregate(data_subset$lma, by = list(data_subset$deciduousness), mean)

# Summary of the model
summary(full_model)
coef_summary <- summary(full_model)$coefficients
saveRDS(full_model, file = "models/full_model_07b_lp.rds")
write.csv(coef_summary, file = "tables/model_coefficients_07b_lp.csv", row.names = TRUE)

## DATA DREDGE REGRESSION
# Load MuMIn package for model selection
library(MuMIn)

# Perform model selection
model_set <- dredge(full_model, trace = TRUE, rank = "AICc")
saveRDS(model_set, file = "models/model_set_07b_lp.rds")
model_set <- readRDS(file = "models/model_set_07b_lp.rds")

# View the top models
head(model_set)

## MODEL AVG REGRESSION
# Model averaging of top models within 2 AICc units
avg_model <- model.avg(model_set, subset = delta < 2)

# Summary of the averaged model
summary(avg_model)
model_summary_output <- capture.output(summary(avg_model))
writeLines(model_summary_output, "tables/model_summary_output_07b_lp.csv")

coef_summary <- summary(avg_model)$coefficients
write.csv(avg_model$msTable, file = "tables/model_fits_07b_lp.csv", row.names = TRUE)
write.csv(coef_summary, file = "tables/modelavg_coefficients_07b_lp.csv", row.names = TRUE)

## VALIDATION
residuals_phylolm <- residuals(full_model)
fitted_values_phylolm <- fitted(full_model)
plot(fitted_values_phylolm, residuals_phylolm,
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

