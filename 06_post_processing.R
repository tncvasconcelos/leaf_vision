# Post-processing results from LM2
# rm(list=ls())
# setwd("~/leaf_computer_vision")
library(data.table)
library(raster)
library(geodata)

#---------------------------------------
# Loading Functions we will need
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
petiole_measurements <- fread("data/width_data_all_specimens_LM2.csv") # lastest results as of Oct 30th 2024
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

# Transform petiole width and leaf area into LMA
merged_dataset$LMA <- NA

# transforming area into m^2 to be easier to compare with other studies
merged_dataset$area <- merged_dataset$area * 0.0001

# transforming petiole width into m as well
merged_dataset$petiole_width <- merged_dataset$petiole_width * 0.01

for(i in 1:nrow(merged_dataset)) {
  one_subset <- subset(merged_dataset, merged_dataset$filename %in% merged_dataset$filename[i])
  petiole_width_for_lma <- min(one_subset$petiole_width)
  merged_dataset$LMA[i] <- make_LMA(merged_dataset$area[i], petiole_width_for_lma) * 100
  cat(i, "of", nrow(merged_dataset), "\r")
}


#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset2.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset2.csv")
#---------------------------------------

# check for per species outliers
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
to_rm_index <- grep("to_rm", colnames(merged_dataset))
merged_dataset <- merged_dataset[rowSums(merged_dataset[,c(66,67,68)]) == 0,]

#---------------------------------------
# Adding latitude and longitude from GBIF metadata
metadata <- list.files(path = "virtual_herbarium_NPleafPaper/metadata/", full.names = T)
metadata <- lapply(metadata, read.csv)
for(i in 1:length(metadata)) {
  if(!"institutionCode" %in% colnames(metadata[[i]])) {
    metadata[[i]]$institutionCode <- NA
  }
  metadata[[i]] <- metadata[[i]][,c("gbifID","decimalLatitude","decimalLongitude","institutionCode")]
  cat(i, "\r")
}
metadata <- do.call(rbind, metadata)

merged_dataset$lat <- NA
merged_dataset$lon <- NA
merged_dataset$inst_code <- NA
for(j in 1:nrow(metadata)) {
  n_lma_result <- grep(metadata$gbifID[j], merged_dataset$filename)
  merged_dataset$lat[n_lma_result] <- metadata$decimalLatitude[j]
  merged_dataset$lon[n_lma_result] <- metadata$decimalLongitude[j]
  merged_dataset$inst_code[n_lma_result] <- metadata$institutionCode[j]
  cat(j, "\r")
}

#write.csv(table(merged_dataset$inst_code), file="inst_codes.csv")

#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset2.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset2.csv")
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

#plot(srad)

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
# write.csv(merged_dataset_test, file="data/merged_dataset.csv", row.names=F)
# merged_dataset_original <- fread("data/old/merged_dataset.csv")
#---------------------------------------

#--------------------------------- 
# By "super-biome":
temperate <- c("Temperate Broadleaf & Mixed Forests", 
               "Temperate Conifer Forests", 
               "Temperate Grasslands, Savannas & Shrublands",
               "Boreal Forests/Taiga",
               "Montane Grasslands & Shrublands",
               "Tundra")
tropical <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", 
              "Tropical & Subtropical Moist Broadleaf Forests", 
              "Tropical & Subtropical Dry Broadleaf Forests", 
              "Tropical & Subtropical Coniferous Forests",
              "Mangroves")
arid <- c("Mediterranean Forests, Woodlands & Scrub",
          "Deserts & Xeric Shrublands")
merged_dataset$super_biome <- NA
for(i in 1:nrow(merged_dataset)) {
  if(merged_dataset$biome[i] %in% temperate) {
    merged_dataset$super_biome[i] <- "temperate"
  } 
  if(merged_dataset$biome[i] %in% tropical) {
    merged_dataset$super_biome[i] <- "tropical"
  }
  if(merged_dataset$biome[i] %in% arid) {
    merged_dataset$super_biome[i] <- "arid"
  }
  cat(i, "\r")
}

############
# By "likely to be decidious":
likely_deciduous <- c("Temperate Broadleaf & Mixed Forests", 
               "Temperate Conifer Forests", 
               "Temperate Grasslands, Savannas & Shrublands",
               "Boreal Forests/Taiga",
               "Tundra",
               "Tropical & Subtropical Dry Broadleaf Forests")
likely_evergreen <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", 
              "Tropical & Subtropical Moist Broadleaf Forests", 
              "Tropical & Subtropical Coniferous Forests",
              "Mangroves","Mediterranean Forests, Woodlands & Scrub",
              "Deserts & Xeric Shrublands",
              "Montane Grasslands & Shrublands")
merged_dataset$deciduousness <- NA
for(i in 1:nrow(merged_dataset)) {
  if(merged_dataset$biome[i] %in% likely_deciduous) {
    merged_dataset$deciduousness[i] <- "likely_deciduous"
  } 
  if(merged_dataset$biome[i] %in% likely_evergreen) {
    merged_dataset$deciduousness[i] <- "likely_evergreen"
  }
  cat(i, "\r")
}


#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset_final.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------

# species_list <- unique(merged_dataset$genus_species)
# results <- c()
# for(i in 1:length(species_list)) {
#   one_subset <- subset(merged_dataset, merged_dataset$genus_species==species_list[i])
#   if(length(unique(one_subset$super_biome))>1) {
#   one_subset$super_biome
#    pp <-  as.data.frame(one_subset[, c("lat","lon")])
#   sp::coordinates(pp) <- ~ lon + lat
#   crs(pp) <- "+proj=longlat +datum=WGS84 +no_defs"
#   points(pp)
#   
#   }
# }



#spot_check <- merged_dataset[,c("genus_species","biome","deciduousness")]
#write.csv(spot_check, file="spot_check_deciduousness.csv", row.names=F)
library(BIEN)
library(data.table)
# Data from BIEN:
list_of_one_trait <- BIEN_trait_trait(trait="whole plant vegetative phenology")
list_of_one_trait0 <- subset(list_of_one_trait, !is.na(list_of_one_trait$scrubbed_species_binomial))
list_of_one_trait0 <- subset(list_of_one_trait0, list_of_one_trait0$access=="public")
list_of_one_trait0 <- list_of_one_trait0[,c("scrubbed_species_binomial","trait_value")]

# Other datasets:
wright_dataset <- as.data.frame(fread("evergreen_deciduous/Wright-etal-2004-LES-supplemental.csv"))
ajb_dataset <- as.data.frame(fread("evergreen_deciduous/ajb216419-sup-0001-appendixs1 (2).csv"))
nph_dataset <- as.data.frame(fread("evergreen_deciduous/nph_3615_sm_tables1-7 (2).csv"))
try_dataset <- as.data.frame(fread("evergreen_deciduous/36766.txt"))

# curating AJB dataset
all_ajb_names <- c()
for(i in 1:length(ajb_dataset$Species)) {
  all_ajb_names <- c(all_ajb_names, paste(unlist(strsplit(ajb_dataset$Species[i]," "))[1:2], collapse=" "))
}
ajb_dataset$Species <- all_ajb_names

# curating TRY dataset
try_dataset <- subset(try_dataset, try_dataset$DataName=="Leaf phenology type")
try_dataset <- try_dataset[,c("AccSpeciesName","OrigValueStr")]
#write.csv(names(table(try_dataset$OrigValueStr)),"evergreen_deciduous/standartizing_try.csv", row.names=F)
standartizing_try <- read.csv("evergreen_deciduous/standartizing_try.csv")
evergreen_ones <- standartizing_try$trait_name[which(standartizing_try$action=="evergreen")] 
deciduous_ones <- standartizing_try$trait_name[which(standartizing_try$action=="deciduous")] 
exclude_ones <- standartizing_try$trait_name[which(standartizing_try$action=="exclude")] 

try_dataset$OrigValueStr[which(try_dataset$OrigValueStr%in%evergreen_ones)] <- "evergreen"
try_dataset$OrigValueStr[which(try_dataset$OrigValueStr%in%deciduous_ones)] <- "deciduous"
try_dataset$OrigValueStr[which(try_dataset$OrigValueStr%in%exclude_ones)] <- "exclude"

try_dataset <- subset(try_dataset, try_dataset$OrigValueStr!="exclude")


# Standartizing column names
colnames(ajb_dataset) <- c("species","leaf_phenology")
colnames(nph_dataset) <- c("species","leaf_phenology")
colnames(wright_dataset) <- c("species","leaf_phenology")
colnames(list_of_one_trait0) <- c("species","leaf_phenology")
colnames(try_dataset) <- c("species","leaf_phenology")

all_pheno_data <- rbind(list_of_one_trait0, ajb_dataset, nph_dataset, wright_dataset, try_dataset)
all_pheno_data <- subset(all_pheno_data, all_pheno_data$leaf_phenology!="")
all_pheno_data <- subset(all_pheno_data, !is.na(all_pheno_data$leaf_phenology))
all_pheno_data <- subset(all_pheno_data, !duplicated(all_pheno_data$species))

all_pheno_data$leaf_phenology[which(all_pheno_data$leaf_phenology=="evergreenergreen")] <- "evergreen"
all_pheno_data <- subset(all_pheno_data, !is.na(all_pheno_data$leaf_phenology))
all_pheno_data <- subset(all_pheno_data, all_pheno_data$leaf_phenology!="variable or conflicting reports")

merged_dataset$genus_species <- gsub("_"," ",merged_dataset$genus_species)
merged_dataset <- merge(as.data.frame(merged_dataset), all_pheno_data, by.x="genus_species",by.y="species",all.x=T)

#max(merged_dataset$LMA)
#which(is.na(merged_dataset$LMA))

#spot_check <- subset(merged_dataset, !is.na(merged_dataset$leaf_phenology))
#unique(spot_check$genus_species)

#---------------------------------------
# Save point
# write.csv(merged_dataset, file="data/merged_dataset_final.csv", row.names=F)
# merged_dataset <- fread("data/merged_dataset.csv")
#---------------------------------------

#boxplot(merged_dataset$LMA~merged_dataset$biome, las=2, horizontal=T)

#boxplot(spot_check$LMA~spot_check$leaf_phenology)
subset_to_check <- subset(merged_dataset, !duplicated(merged_dataset$genus_species))
write.csv(subset_to_check[,c("genus_species","leaf_phenology")], file="subset_to_check.csv", row.names = F)


# merged_dataset_original <- as.data.frame(merged_dataset_original)[,c(1,70:ncol(merged_dataset_original))]
# merged_dataset_test <- merge(merged_dataset, merged_dataset_original, by="component_name", all=T)
# merged_dataset_test <- subset(merged_dataset_test, !is.na(merged_dataset_test$bio_1))
