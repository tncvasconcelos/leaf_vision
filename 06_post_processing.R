setwd("/Users/tvasc/Desktop/leaf_computer_vision")
setwd("~/leaf_vision/")
#rm(list=ls())
library(data.table)
library(raster)

# WWFload is taken from speciesgeocodeR; all credit goes to the original authors
WWFload <- function(x = NULL) {
  if (missing(x)) {
    x <- getwd()
  }
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
  unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
  file.remove(file.path(x, "wwf_ecoregions.zip"))
  wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
                                              "wwf_terr_ecos.shp"))
  return(wwf)
}

localityToBiome <- function (points, lat="lat",lon="lon") {
  #colnames(points) <- c("acceptedScientificName","key","decimalLatitude","decimalLongitude","basisOfRecord","issues")
  cat("Getting biome from locality data...")
  points[,lat] <-  as.numeric(points[,lat])
  points[,lon] <-  as.numeric(points[,lon])
  locations.spatial <- sp::SpatialPointsDataFrame(coords=points[,c(which(colnames(points)==lon), which(colnames(points)==lat))], data=points)
  wwf <- WWFload(tempdir())
  mappedregions <- sp::over(locations.spatial, wwf)
  biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
  points$eco_name <- mappedregions$ECO_NAME
  points$biome <- biomes[mappedregions$BIOME]
  return(points)
}

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
  cat(i, "of", nrow(merged_dataset), "\r")
}

write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)

#-------------------
# merged_dataset <- fread("data/merged_dataset.csv")

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
  cat(i, "of", nrow(merged_dataset), "\r")
}

write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)

#-------------------
# Create a mean LMA for each specimen 
#merged_dataset <- fread("data/merged_dataset.csv")
spp <- unique(merged_dataset$genus_species)
merged_dataset$to_rm_lma <- 0
merged_dataset$to_rm_wdth <- 0
merged_dataset$to_rm_area <- 0
for(i in seq_along(spp)){
  cat("\r", i, "of", length(spp))
  tmp <- merged_dataset[merged_dataset$genus_species == spp[i],]
  to_rm <- tmp$component_id[abs((mean(log(tmp$LMA)) - 
      log(tmp$LMA))/sd(log(tmp$LMA))) > 3]
  if(length(to_rm) >= 1 & nrow(tmp) > 1){
    merged_dataset$to_rm_lma[match(to_rm, merged_dataset$component_id)] <- 1
  }
  to_rm <- tmp$component_id[abs((mean(log(tmp$petiole_width)) - 
      log(tmp$petiole_width))/sd(log(tmp$petiole_width))) > 3]
  if(length(to_rm) >= 1 & nrow(tmp) > 1){
    merged_dataset$to_rm_wdth[match(to_rm, merged_dataset$component_id)] <- 1
  }
  to_rm <- tmp$component_id[abs((mean(log(tmp$area)) - 
      log(tmp$area))/sd(log(tmp$area))) > 3]
  if(length(to_rm) >= 1 & nrow(tmp) > 1){
    merged_dataset$to_rm_area[match(to_rm, merged_dataset$component_id)] <- 1
  }
}

merged_dataset <- merged_dataset[rowSums(merged_dataset[,65:67]) == 0,]

#-------------------
# Adding lat and lon
metadata <- list.files(path = "virtual_herbarium_NPleafPaper/metadata/", full.names = T)
metadata <- lapply(metadata, read.csv)
for(i in 1:length(metadata)) {
  metadata[[i]] <- metadata[[i]][,c("gbifID","decimalLatitude","decimalLongitude")]
  cat(i, "\r")
}
metadata <- do.call(rbind, metadata)

merged_dataset$lat <- NA
merged_dataset$lon <- NA
for(j in 1:nrow(metadata)) {
  n_lma_result <- grep(metadata$gbifID[j], merged_dataset$filename)
  merged_dataset$lat[n_lma_result] <- metadata$decimalLatitude[j]
  merged_dataset$lon[n_lma_result] <- metadata$decimalLongitude[j]
  cat(j, "\r")
}
write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)

#-------------------
# Adding climate data
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$lon))
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$lat))

#---------------------------------
# Adding biome data
biomes_for_points <- localityToBiome(as.data.frame(merged_dataset), lat="lat",lon="lon")
merged_dataset <- cbind(merged_dataset, biomes_for_points)

#---------------------------------
# Extracting climate and altitude medians per ecoregion
bio <- raster::getData('worldclim', var='bio', res=2.5) # 19 worldclim vars
#alt <- raster::getData('worldclim', var='alt', res=2.5) # altitude
ai <- raster("layers/ai_v3_yr.tif")
npp <- raster("layers/MOD17A3H_Y_NPP_2023-01-01_rgb_720x360.TIFF")
wind <- raster("layers/mean_wind.tif")
srad <- raster("layers/mean_srad.tif")
et0 <- raster("layers/et0_v3_yr.tif")

#---------------------------------
coordinates <- merged_dataset[,c("lat","lon")]

points <- coordinates
layers <- c(as.list(bio), ai, npp, wind, srad, et0)
names(layers) <- c(names(bio),"ai","npp","wind","srad","et0")

for(layer_index in 13:length(layers)) {
  layer <- layers[[layer_index]]
  all_values <- c()
  for(i in 1:nrow(points)) {
    one_point <- points[i,]
    one_point <- as.data.frame(one_point)
    sp::coordinates(one_point) <- ~ lon + lat
    crs(one_point) <- "+proj=longlat +datum=WGS84 +no_defs"
    #buffered_point <- buffer(one_point, width = 10000)
    values <- raster::extract(layer, one_point)
    all_values <- c(all_values, values)
    cat(i, "\r")
  }  
  coordinates <- cbind(coordinates, all_values)
  colnames(coordinates)[2+layer_index] <- names(layers)[layer_index]
}
merged_dataset <- cbind(merged_dataset,coordinates)

# removing 0 coordinates because I didnt do this before
merged_dataset <- subset(merged_dataset, merged_dataset$lat!=0 & merged_dataset$lat!=0)
write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)

# plot(log(merged_dataset$area)~log(merged_dataset$bio12))
# test <- lm(log(merged_dataset$area)~log(merged_dataset$bio12))
# summary(test)
# abline(test, col="red")

lma_results <- aggregate(merged_dataset$LMA, list(merged_dataset$genus_species), 
  FUN = function(x) c(mean(log(x)), sd(log(x))/length(x)))
lma_results <- data.frame(sp = lma_results$Group.1,
  lma = lma_results$x[,1],
  se =  lma_results$x[,2])

# 
# for(i in 1:nrow(lma_results)) {
#   one_specimen <- subset(merged_dataset, merged_dataset$filename==lma_results$specimen[i])
#   lma_results$mean_LMA[i] <- mean(one_specimen$LMA)
#   cat(i, "\r")
# }

pdf("results/LMA_dist.pdf")
hist(lma_results$lma, breaks=100, xlab="log(LMA)", main="LMA distribution")
dev.off()
write.csv(lma_results, file="data/lma_results.csv", row.names=F)

# merged_dataset <- subset(merged_dataset, !merged_dataset$petiole_width %in% tail(sort(merged_dataset$petiole_width), n=1000))
# merged_dataset <- subset(merged_dataset, !merged_dataset$petiole_width %in% head(sort(merged_dataset$petiole_width), n=1000))
# 
# quantile(merged_dataset$petiole_width, probs = 0.75)

#plot(log(lma_results$mean_LMA)~lma_results$lat)
write.csv(lma_results, file="lma_results.csv", row.names=F)


