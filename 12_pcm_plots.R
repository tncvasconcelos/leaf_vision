#------------------------------------

# Load necessary libraries
library(phylolm)
library(ggplot2)
library(dplyr)

load("results/data_subset_w_leaf_phenology.Rsave")

model <- phylolm(lma~bio_1, data=dat, phy=phy)
summary(model)

plot(lma~bio_1,data=dat)
abline(model, col="red")

# Generate a sequence of bio_5 values for smooth lines
bio_5_seq <- seq(min(data_subset$bio_5), max(data_subset$bio_5), length.out = 100)

# Create a data frame for each deciduousness level for predictions
pred_evergreen <- data.frame(bio_5 = bio_5_seq, deciduousness = "evergreen")
pred_deciduous <- data.frame(bio_5 = bio_5_seq, deciduousness = "deciduous")

# Predict lma for each level using the phylogenetic model
pred_evergreen$lma <- predict(model_evergreen, newdata = pred_evergreen)
pred_deciduous$lma <- predict(model_deciduous, newdata = pred_deciduous)

# Combine the prediction data
pred_data <- rbind(pred_evergreen, pred_deciduous)


# Plot with solid lines
bio_5_plot <- ggplot(data_subset, aes(x = bio_5, y = lma, color = deciduousness)) +
  geom_point(size = 1) +
  geom_line(data = pred_data, aes(y = lma, color = deciduousness), size = 0.8) +
  labs(color = "Deciduousness") +
  theme_bw() +
  ggtitle("Phylogenetically Corrected Regression of LMA on BIO_5") +
  scale_color_manual(values = c("evergreen" = "darkblue", "deciduous" = "orange"))

#----------------------------------


# Fit phylogenetic regression models for each level of deciduousness
model_evergreen <- phylolm(lma ~ ai, data = filter(data_subset, deciduousness == "evergreen"), phy = phy)
model_deciduous <- phylolm(lma ~ ai, data = filter(data_subset, deciduousness == "deciduous"), phy = phy)

# Create predictions for plotting
data_subset <- data_subset %>%
  mutate(pred_evergreen = ifelse(deciduousness == "evergreen", predict(model_evergreen, newdata = .), NA),
         pred_deciduous = ifelse(deciduousness == "deciduous", predict(model_deciduous, newdata = .), NA))

# Generate a sequence of bio_5 values for smooth lines
ai_seq <- seq(min(data_subset$ai), max(data_subset$ai), length.out = 100)

# Create a data frame for each deciduousness level for predictions
pred_evergreen <- data.frame(ai = ai_seq, deciduousness = "evergreen")
pred_deciduous <- data.frame(ai = ai_seq, deciduousness = "deciduous")

# Predict lma for each level using the phylogenetic model
pred_evergreen$lma <- predict(model_evergreen, newdata = pred_evergreen)
pred_deciduous$lma <- predict(model_deciduous, newdata = pred_deciduous)

# Combine the prediction data
pred_data <- rbind(pred_evergreen, pred_deciduous)

# Plot with solid lines
bio_5_plot <- ggplot(data_subset, aes(x = ai, y = lma, color = deciduousness)) +
  geom_point(size = 1) +
  geom_line(data = pred_data, aes(y = lma, color = deciduousness), size = 0.8) +
  labs(color = "Deciduousness") +
  theme_bw() +
  ggtitle("Phylogenetically Corrected Regression of LMA on AI") +
  scale_color_manual(values = c("evergreen" = "darkblue", "deciduous" = "orange"))



