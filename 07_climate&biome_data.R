setwd("/Users/tvasc/Desktop/leaf_computer_vision")
# rm(list=ls())
library(data.table)

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

lma_results <- read.csv("lma_results.csv")
lma_results <- subset(lma_results, !is.na(lma_results$lon))

#---------------------------------
# Adding biome data
biomes_for_points <- localityToBiome(lma_results, lat="lat",lon="lon")
lma_results <- cbind(lma_results, biomes_for_points)

#lma_results <- subset(lma_results, !lma_results$mean_LMA %in% boxplot_lma$out)
#boxplot(lma_results$mean_LMA~lma_results$biome,las = 2)

#---------------------------------
# Binarizing biome type:
#closed_canopy <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga")
#open_canopy <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra","Deserts & Xeric Shrublands", "Mediterranean Forests, Woodlands & Scrub",  "Mangroves")
#lma_results$bin_biome <- unlist(ifelse(lma_results$biome%in%closed_canopy, "closed","open"))


#boxplot(lma_results$mean_LMA~lma_results$bin_biome)


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
coordinates <- lma_results[,c("lat","lon")]

points <- coordinates
layers <- c(as.list(bio), ai, npp, wind, srad, et0)
names(layers) <- c(names(bio),"ai","npp","wind","srad","et0")

for(layer_index in 1:length(layers)) {
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
lma_results <- cbind(lma_results,coordinates[,3:4])

plot(log(lma_results$mean_LMA)~log(lma_results$bio2))

#--------------------------------- 
# By "super-biome":
temperate <- c("Temperate Broadleaf & Mixed Forests", 
               "Temperate Conifer Forests", 
               "Temperate Grasslands, Savannas & Shrublands",
               "Boreal Forests/Taiga",
               "Montane Grasslands & Shrublands", 
               "Tropical & Subtropical Coniferous Forests",
               "Tundra")
tropical <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", 
              "Tropical & Subtropical Moist Broadleaf Forests", 
              "Tropical & Subtropical Dry Broadleaf Forests", 
              "Mangroves")
arid <- c("Mediterranean Forests, Woodlands & Scrub",
          "Deserts & Xeric Shrublands")
all_surveys$super_biome <- NA
for(i in 1:nrow(all_surveys)) {
  if(all_surveys$biome[i] %in% temperate) {
    all_surveys$super_biome[i] <- "temperate"
  } 
  if(all_surveys$biome[i] %in% tropical) {
    all_surveys$super_biome[i] <- "tropical"
  }
  if(all_surveys$biome[i] %in% arid) {
    all_surveys$super_biome[i] <- "arid"
  }
}

#--------------------------------- 
# Binarizing tropical vs. temperate:
tropical <- which(all_surveys$latitude > -23 & all_surveys$latitude < 23) 
temperate <- which(all_surveys$latitude <= -23 | all_surveys$latitude >= 23) 
all_surveys$tropical <- NA
all_surveys$tropical[tropical] <- "tropical"
all_surveys$tropical[temperate] <- "temperate"


#--------------------------------- 
saveRDS(all_surveys, "data/community_studies_w_habitat_categories_&_env_vars.Rdata")
write.csv(all_surveys, "data/community_studies_w_habitat_categories_&_env_vars.csv")