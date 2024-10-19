setwd("/Users/tvasc/Desktop/leaf_computer_vision")
#rm(list=ls())
library(data.table)

# function to make LMA from leaf area and petiole width
# for woody dicots (based on Royer et al. 2007, r^2=0.55)
make_LMA <- function(leaf_area, petiole_width) {
  LMA = 3.07 + 0.382 * log(petiole_width^2 / leaf_area)
  return(exp(LMA))
}

# Read Leaf Machine output
all_measurements <- fread("data/LM2_MEASUREMENTS_CLEAN.csv")
# Read petiole width data
petiole_measurements <- fread("data/width_data.csv")
# Merge datasets
merged_dataset <- merge(all_measurements, petiole_measurements, by.x="component_name",by.y="filename")
# remove NAs from area data
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$area))

merged_dataset$petiole_width <- NA
for(i in 1:nrow(merged_dataset)) {
  one_petiole <- merged_dataset$width_pixels[i]
  if(merged_dataset$conversion_mean[i]!=0){
    merged_dataset$petiole_width[i] <- one_petiole / merged_dataset$conversion_mean[i]
  } else {
    merged_dataset$petiole_width[i] <- one_petiole / merged_dataset$predicted_conversion_factor_cm[i]
  }
  cat(i, "\r")
}

write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)

###################
###################
###################
merged_dataset <- fread("data/merged_dataset.csv")

# some summary stats:
min(merged_dataset$area)
merged_dataset$filename[which.min(merged_dataset$area)]
max(merged_dataset$area)
merged_dataset$filename[which.max(merged_dataset$area)]

pdf("results/leaf_area_dist.pdf")
hist(log(merged_dataset$area), breaks=100, xlab="log(leaf area cm^2)", main="Leaf area distribution")
dev.off()

min(merged_dataset$petiole_width)
merged_dataset$filename[which.min(merged_dataset$petiole_width)]
max(merged_dataset$petiole_width)
merged_dataset$filename[which.max(merged_dataset$petiole_width)]

pdf("results/petiole_width_dist.pdf")
hist(log(merged_dataset$petiole_width), breaks=100, xlab="log(petiole width cm)", main="Petiole width distribution")
dev.off()

pdf("results/cor_leaf_area_petiole_width.pdf")
model <- lm(log(merged_dataset$area)~log(merged_dataset$petiole_width))
summary(model)
plot(log(merged_dataset$area)~log(merged_dataset$petiole_width), 
     xlab="log(petiole width cm)", ylab="log(leaf area cm^2)")
abline(model, col="red")
dev.off()


# Transform numbers into LMA
merged_dataset$LMA <- NA
for(i in 1:nrow(merged_dataset)) {
  merged_dataset$LMA[i] <- make_LMA(merged_dataset$area[i], merged_dataset$petiole_width[i])
  cat(i, "\r")
}

write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)

# Create a mean LMA for each specimen 

lma_results <- data.frame(specimen=unique(merged_dataset$filename), mean_LMA=NA)
for(i in 1:nrow(lma_results)) {
  one_specimen <- subset(merged_dataset, merged_dataset$filename==lma_results$specimen[i])
  lma_results$mean_LMA[i] <- mean(one_specimen$LMA)
  cat(i, "\r")
}

pdf("results/LMA_dist.pdf")
hist(log(lma_results$mean_LMA), breaks=100, xlab="log(LMA)", main="LMA distribution")
dev.off()
write.csv(lma_results, file="lma_results.csv", row.names=F)

# merged_dataset <- subset(merged_dataset, !merged_dataset$petiole_width %in% tail(sort(merged_dataset$petiole_width), n=1000))
# merged_dataset <- subset(merged_dataset, !merged_dataset$petiole_width %in% head(sort(merged_dataset$petiole_width), n=1000))
# 
# quantile(merged_dataset$petiole_width, probs = 0.75)

metadata <- list.files(path = "virtual_herbarium_NPleafPaper/metadata/", full.names = T)
metadata <- lapply(metadata, read.csv)
for(i in 1:length(metadata)) {
  metadata[[i]] <- metadata[[i]][,c("gbifID","decimalLatitude","decimalLongitude")]
  cat(i, "\r")
}
metadata <- do.call(rbind, metadata)

lma_results$lat <- NA
lma_results$lon <- NA
for(j in 1:nrow(metadata)) {
  n_lma_result <- grep(metadata$gbifID[j], lma_results$specimen)
  lma_results$lat[n_lma_result] <- metadata$decimalLatitude[j]
  lma_results$lon[n_lma_result] <- metadata$decimalLongitude[j]
  cat(j, "\r")
}

#plot(log(lma_results$mean_LMA)~lma_results$lat)
write.csv(lma_results, file="lma_results.csv", row.names=F)


