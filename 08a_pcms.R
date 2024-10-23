#unstandardized values of LA
library(ape)
library(phylolm)
library(geiger)
library(phytools)
library(scales)
library(RColorBrewer)

setwd("~/leaf_vision/")

merged_dataset <- read.csv("data/merged_dataset.csv")
tre <- read.tree("trees/GBMB.tre")

la_results <- aggregate(merged_dataset$area, list(merged_dataset$genus_species), 
  FUN = function(x) c(mean(log(x)), sd(log(x))/length(x)))
la_results <- data.frame(sp = la_results$Group.1,
  la = la_results$x[,1],
  se =  la_results$x[,2])

merged_dataset_2 <- merged_dataset[!duplicated(merged_dataset$filename),]
focal_cols <- grep("bio_1$", colnames(merged_dataset_2)):ncol(merged_dataset_2)
climate_data <- aggregate(merged_dataset_2[,focal_cols], list(merged_dataset_2$genus_species), 
  FUN = function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm = TRUE)/length(na.omit(x))))
climate_data <- do.call(cbind, lapply(climate_data, function(x) as.data.frame(x)))[,-1]
colnames(climate_data) <- gsub("\\.V1", "", colnames(climate_data))
colnames(climate_data) <- gsub("\\.V2", "_se", colnames(climate_data))
dat <- cbind(la_results, climate_data)

missing_sp <- dat$sp[!dat$sp %in% tre$tip.label]
dat <- dat[dat$sp %in% tre$tip.label, ]

phy <- keep.tip(tre, dat$sp)
phy <- ladderize(phy)
phy$node.label <- NULL
rownames(dat) <- dat$sp
dat <- dat[phy$tip.label,]
H <- max(node.depth.edgelength(phy))

write.csv(dat, file = "data/la_dat_a.csv", row.names = FALSE)

rescale_values <- function(x, a, b) {
  (a + (x - min(x)) / (max(x) - min(x)) * (b - a))
}
plot_tip_values <- setNames(rescale_values(dat$la, H*1.01, H*1.1), dat$sp)

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
ggplot(dat, aes(x = la)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Leaf Area (LA)", x = "Leaf Area", y = "Frequency")

# Melt the data
library(reshape2)
predictor_vars <- grep("_se$", names(dat), invert = TRUE, value = TRUE)
predictor_vars <- predictor_vars[!predictor_vars %in% c("sp", "la", "se")]
predictor_vars <- predictor_vars[predictor_vars!="et0"]


## CORRELATION EXAMINATION
# Compute correlation matrix
cor_mat <- cor(dat[, predictor_vars], use = "complete.obs")

library(caret)
library(phylolm)
# Find variables with correlation higher than 0.7
high_cor <- findCorrelation(cor_mat, cutoff = 0.7)
reduced_vars <- predictor_vars[-high_cor]

# Updated list of variables
print(reduced_vars)

## PCA EXAMINATION
# Perform PCA on standardized variables
pca_res <- prcomp(dat[, predictor_vars], scale. = TRUE)

# Scree plot to determine the number of principal components to retain
plot(pca_res, type = "l")

# Get the proportion of variance explained
cumsum(summary(pca_res)$importance[2, ])

## FULL PHYLO REGRESSION
# Create a comparative data object
# comp_data <- comparative.data(phy, dat[,c("sp", "la", reduced_vars)], 
#   names.col = "sp", vcv = TRUE, na.omit = TRUE)

# Full model with all predictors (after variable selection)
# full_model <- pgls(lma ~ ., data = comp_data, lambda = "ML")
data_subset <- dat[, c("la", reduced_vars)]
formula_full <- as.formula(paste("la ~", paste(reduced_vars, collapse = " + ")))
full_model <- phylolm(formula_full, phy = phy, data = dat, model = "lambda", REML = FALSE)

# Summary of the model
summary(full_model)

## DATA DREDGE REGRESSION
# Load MuMIn package for model selection
library(MuMIn)

# Perform model selection
model_set <- dredge(full_model, trace = TRUE, rank = "AICc")

# View the top models
head(model_set)

## MODEL AVG REGRESSION
# Model averaging of top models within 2 AICc units
avg_model <- model.avg(model_set, subset = delta < 2)

# Summary of the averaged model
summary(avg_model)

## VALIDATION
residuals_phylolm <- residuals(full_model)
fitted_values_phylolm <- fitted(full_model)
plot(fitted_values_phylolm, residuals_phylolm,
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")


