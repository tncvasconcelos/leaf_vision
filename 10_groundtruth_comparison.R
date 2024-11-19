# setwd("/Users/tvasc/Desktop/leaf_computer_vision")
# setwd("~/leaf_vision/")
# rm(list=ls())
library(data.table)

# load ground truthing data
manual_measurements <- as.data.frame(fread("data/GT_comparison.csv"))

pdf("plots/gt_comparison.pdf")
plot(manual_measurements$width_pixels~manual_measurements$pixel_distance, 
     xlim=c(0,85), ylim=c(0,85),
     xlab="petiole width in pixels (LM2)",
     ylab="petiole width in pixels (manual)")
abline(a=0,b=1, col="red")
dev.off()

comparison <- lm(manual_measurements$width_pixels~manual_measurements$pixel_distance)
summary(comparison)


