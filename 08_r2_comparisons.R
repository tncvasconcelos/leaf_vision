library(ape)
library(MuMIn)
library(rr2)
library(phylolm)
setwd("~/leaf_vision/")

# 07a
full_model <- readRDS(file = "models/full_model_07b.rds")
model_set <- readRDS(file = "models/model_set_07b.rds")
phy <- read.tree("trees/phy_070a.tre")
data_subset <- dat <- full_model$model.frame
full_formula <- full_model$formula

##### fitting an lm
lm_model <- lm(full_formula, data = dat, na.action = na.fail)
model_set_lm <- dredge(lm_model, trace = TRUE, rank = "AICc")
model_avg_lm <- model.avg(model_set_lm, subset = delta < 2)
avg_coefs <- coef(model_avg_lm)
formula_avg <- as.formula(paste("lma ~", paste(names(avg_coefs)[-1], collapse = " + ")))
avg_model_lm <- lm(formula_avg, data = dat)

##### Rsqared based on refitting model average
# Perform model averaging with MuMIn
model_avg <- model.avg(model_set, subset = delta < 2)
# Extract averaged coefficients
avg_coefs <- coef(model_avg)
# Create a formula for the new model with averaged coefficients
formula_avg <- as.formula(paste("lma ~", paste(names(avg_coefs)[-1], collapse = " + ")))
# Run a new phylolm model using the averaged coefficients
avg_model <- phylolm(formula_avg, data = dat, phy = phy, model = "lambda", REML = FALSE)
# Summary to extract R-squared
summary_avg_model <- summary(avg_model)
r_squared_avg <- summary_avg_model$r.squared
r_squared_avg

##### Rsqared based on model averaging results
# Extract R-squared values and weights for each model
index <- which(model_set$delta < 2)
r_squared_values <- sapply(attr(model_set[index], "model.calls"), function(mod_call) {
  mod <- eval(mod_call)
  mod$r.squared
})
model_weights <- model_set$weight[index]/sum(model_set$weight[index])
# Calculate weighted R-squared
weighted_r_squared <- sum(r_squared_values * model_weights, na.rm = TRUE)
weighted_r_squared

#### IVES rr2
avg_model_phy <- phylolm(formula_avg, data = dat, phy = phy, model = "lambda", REML = FALSE)
R2_resid(avg_model_phy, phy = phy)
R2_lik(avg_model_phy)
R2_resid(avg_model_lm)
R2_lik(avg_model_lm)

