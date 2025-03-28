# rm(list=ls())
# setwd("~/Desktop/leaf_computer_vision")
# setwd("~/leaf_vision/")
library(gridExtra)
library(data.table)
library(ggplot2)

#---------------------------------------
merged_dataset <- read.csv("data/merged_dataset_final.csv")
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$area))
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$petiole_width))

# gt petiole comparison
manual_measurements <- as.data.frame(fread("data/GT_comparison.csv"))
model <- lm(manual_measurements$width_pixels~manual_measurements$pixel_distance)

# Acer

manual_measurements$genus <- unlist(lapply(strsplit(manual_measurements$filename,"_"),"[",1))

all_genera <- unique(manual_measurements$genus)
r_squared <- c()
for(i in 1:length(all_genera)) {
  one_genus_subset <- subset(manual_measurements, manual_measurements$genus==all_genera[i])
  if(nrow(one_genus_subset)>10) {
    model <- lm(one_genus_subset$width_pixels~ one_genus_subset$pixel_distance)
    r_squared <- c(r_squared, summary(model)$r.squared)
    names(r_squared)[length(r_squared)] <- all_genera[i]
  }
}
sort(r_squared, decreasing=T)

write.csv(data.frame(r_squared), file="lm_ten_specimens_error.csv")

