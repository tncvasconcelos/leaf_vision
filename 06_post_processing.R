# setwd("/Users/tvasc/Desktop/leaf_computer_vision")
# setwd("~/leaf_vision/")
# rm(list=ls())
library(data.table)
library(raster)
library(geodata)

#---------------------------------------
# Functions
#---------------------------------------
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

# Get biome from coordinate
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

# Function to make LMA from leaf area and petiole width
# for woody dicots (based on Royer et al. 2007, r^2=0.55)
make_LMA <- function(leaf_area, petiole_width) {
  LMA = 3.07 + 0.382 * log(petiole_width^2 / leaf_area)
  return(exp(LMA))
}

# Function to filter wrong coordinates based on POWO
FilterWCVP_genus <- function(points, all_vars, twgd_data, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  dubious_points <- c()
  all_genera <- unique(all_vars$genus)
  all_vars_genus_level <- all_vars[all_vars$taxon_rank=="Genus",]
  for(genus_index in 1:length(all_genera)) {
    gbif_subset <- subset(tmp_points, tmp_points$genus == all_genera[genus_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars_genus_level, all_vars_genus_level$genus == all_genera[genus_index])
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      if(nrow(area_plus_buffer)>0) {
        coords <- gbif_subset[,c("x","y")]
        sp::coordinates(coords) <- ~ x + y
        answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
        if(length(answer) != 0) {
          dubious_points <- c(dubious_points, as.character(gbif_subset$filename[answer]))
        }
      }
    }
    cat(genus_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$filename) %in% dubious_points)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
} 
#---------------------------------------
#---------------------------------------
#---------------------------------------

# Read Leaf Machine output
all_measurements <- fread("data/LM2_MEASUREMENTS_CLEAN.csv") #this is on gitignore due to size
# Read petiole width data
petiole_measurements <- fread("data/width_data.csv")
# Merge datasets
merged_dataset <- merge(all_measurements, petiole_measurements, by.x="component_name",by.y="filename")
# remove NAs from area data
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$area))

# Transform petiole width from pixels into metric system using conversion factors
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
#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------

# Some summary stats and plots:
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

#---------------------------------------
# Transform petiole width and leaf area into LMA
merged_dataset$LMA <- NA
for(i in 1:nrow(merged_dataset)) {
  merged_dataset$LMA[i] <- make_LMA(merged_dataset$area[i], merged_dataset$petiole_width[i])
  cat(i, "of", nrow(merged_dataset), "\r")
}

#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------

# Create a mean LMA for each specimen 
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

#---------------------------------------
# Adding latitude and longite from GBIF metadata
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

#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------
# Filtering wrong coordinates from POWO (also should've done this before)
dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# Looking at the POWO table and TDWG to clean occurrence points
path="wcvp/wgsrpd-master/level3/level3.shp"
twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))

merged_dataset <- FilterWCVP_genus(points=merged_dataset, all_vars, twgd_data, lon="lon",lat="lat")

#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------

#---------------------------------------
# Adding climate and biome data
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$lon))
merged_dataset <- subset(merged_dataset, !is.na(merged_dataset$lat))
merged_dataset <- localityToBiome(as.data.frame(merged_dataset), lat="lat",lon="lon")

#---------------------------------------
# Extracting climate, altitude and other things for each points
bio <- worldclim_global(var="bio", res=2.5, version="2.1", path=getwd())# 19 worldclim vars
names(bio) <- gsub("wc2.1_2.5m_","",names(bio))
bio_clim <- list()
for(i in 1:length(names(bio))){
  bio_clim[[i]] <- raster(bio[[i]])
  names(bio_clim)[i] <- names(bio)[i]
}
alt <- raster(elevation_global(res=2.5, path=getwd()))# altitude
wind <- raster("layers/mean_wind.tif")
srad <- raster("layers/mean_srad.tif")
et0 <- raster("layers/et0_v3_yr.tif")
ai <- raster("layers/ai_v3_yr.tif")

coordinates <- merged_dataset[,c("lat","lon")]
points <- coordinates
layers <- c(bio_clim, alt, ai, wind, srad, et0)
names(layers) <- c(names(bio_clim),"alt", "ai","wind","srad","et0")

for(layer_index in 1:length(layers)) {
  layer <- layers[[layer_index]]
  all_values <- c()
  for(i in 1:nrow(points)) {
    one_point <- points[i,]
    one_point <- as.data.frame(one_point)
    sp::coordinates(one_point) <- ~ lon + lat
    crs(one_point) <- "+proj=longlat +datum=WGS84 +no_defs"
    values <- raster::extract(layer, one_point)
    all_values <- c(all_values, values)
    cat(i, "\r")
  }  
  coordinates <- cbind(coordinates, all_values)
  colnames(coordinates)[2+layer_index] <- names(layers)[layer_index]
}

merged_dataset <- cbind(merged_dataset,coordinates[,3:ncol(coordinates)])


#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------

#plot(log(merged_dataset$area)~log(merged_dataset$))
# test <- lm(merged_dataset$LMA~merged_dataset$srad)
# summary(test)
# abline(test, col="red")

# 
# for(i in 1:nrow(lma_results)) {
#   one_specimen <- subset(merged_dataset, merged_dataset$filename==lma_results$specimen[i])
#   lma_results$mean_LMA[i] <- mean(one_specimen$LMA)
#   cat(i, "\r")
# }

# pdf("results/LMA_dist.pdf")
# hist(lma_results$lma, breaks=100, xlab="log(LMA)", main="LMA distribution")
# dev.off()
# write.csv(lma_results, file="data/lma_results.csv", row.names=F)

# merged_dataset <- subset(merged_dataset, !merged_dataset$petiole_width %in% tail(sort(merged_dataset$petiole_width), n=1000))
# merged_dataset <- subset(merged_dataset, !merged_dataset$petiole_width %in% head(sort(merged_dataset$petiole_width), n=1000))
# 
# quantile(merged_dataset$petiole_width, probs = 0.75)

#plot(log(lma_results$mean_LMA)~lma_results$lat)
#write.csv(lma_results, file="lma_results.csv", row.names=F)


