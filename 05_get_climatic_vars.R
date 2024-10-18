
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


# getting biomes for each species
getBiomes <- function (points, species="species") {
  cat("Summarizing biome from locality data...")
  points <- as.data.frame(points) # not sure how to do it without transforming back to data.frame
  points <- subset(points, !is.na(points[,"biome"]))
  categories <- unique(points[,"biome"])
  taxa <- as.character(unique(points[,species]))
  result <- matrix(0, nrow=length(taxa), ncol=length(categories))
  rownames(result) <- taxa
  colnames(result) <- categories
  cat("\n")
  for (taxon_index in seq_along(taxa)) {
    for (category_index in seq_along(categories)) {
      x0 <- points[,species]==taxa[taxon_index]
      x1 <- points[,"biome"]==categories[category_index]
      result[taxon_index, category_index] <- length(which(x0 & x1))
    }
    cat(taxon_index, "\r")
  }
  return(result)
}

setwd("/Users/tvasc/Desktop/leaf_computer_vision")

all_metadata <- list.files("virtual_herbarium_NPleafPaper/metadata/", full.names = T)
all_metadata <- lapply(all_metadata, read.csv)

decimalLatitude <- unlist(lapply(all_metadata, "[[", "decimalLatitude"))
decimalLongitude <- unlist(lapply(all_metadata, "[[", "decimalLongitude"))
gbifID <- unlist(lapply(all_metadata, "[[", "gbifID"))
scientificName <- unlist(lapply(all_metadata, "[[", "scientificName"))

occ_data <- data.frame(scientificName,gbifID,decimalLatitude,decimalLongitude)
occ_data_biomes <- localityToBiome(occ_data, lat="decimalLatitude",lon="decimalLongitude")

sorted_counts <- sort(table(occ_data_biomes$biome), decreasing=T)
par(mar = c(6, 4, 4, 2) + 0.1)
bar_positions <- barplot(sorted_counts, main = "Specimens by Institution",
                         xlab = "", ylab = "Number of Specimens",
                         col = "forestgreen", las = 1, names.arg = NA, cex.main = 1, cex.lab = 0.75)
text(x = bar_positions, y = par("usr")[3] - (par("usr")[4]*.05), labels = names(sorted_counts),
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)
par(mar = c(5, 4, 4, 2) + 0.1)

#head(sort(table(occ_data_biomes$eco_name),decreasing=T))




