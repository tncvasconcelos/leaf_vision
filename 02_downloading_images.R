# rm(list=ls())
library(rgbif)
library(magick)
source("00_functions.R")
setwd("/Users/tvasc/Desktop/leaf_computer_vision")

# If local
path="TWDG/wgsrpd-master/level3/level3.shp"
#-----------------------------
twgd_data <- maptools::readShapeSpatial(path)
twgd_data01 <- sf::st_as_sf(twgd_data)

#-----------------------------
# Load WCVP dataset
dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_vars <- subset(all_vars, all_vars$taxon_rank=="Species")
all_vars <- subset(all_vars, all_vars$taxon_status=="Accepted")

#-----------------------------
all_available_data <- read.csv("available_data.csv")
all_available_data <- subset(all_available_data, !is.na(all_available_data$V2))
all_available_data <- subset(all_available_data, all_available_data$V2 >= 5)


all_species_to_sample <- all_available_data$V1
for(species_index in 1:length(all_species_to_sample)) { 
  metadata <- NULL
  try(metadata <- full.herb.search(all_species_to_sample[species_index]), silent=T)
  if(!is.null(metadata)) {
    #--------------------------------------
    # Filtering dataset little before downloading
    # remove things with no coordinates
    metadata <- subset(metadata, !is.na(metadata$decimalLatitude))
    metadata <- subset(metadata, !is.na(metadata$decimalLongitude))
    
    # removing inaturalist observations
    metadata <- subset(metadata, metadata$basisOfRecord=="PRESERVED_SPECIMEN")
    metadata <- subset(metadata, !grepl("inaturalist",metadata$media_url))
    
    #--------------------------------------        
    # removing records out of the native range for the genus
    #--------------------------------------    
    metadata <- FilterWCVP_genus(points=metadata, all_vars, twgd_data, lon="decimalLongitude", lat="decimalLatitude") 

    # if more than 10 observation by the end of the filtering keep only 10
    if(nrow(metadata) >= 10) {
      metadata <- metadata[sample(1:nrow(metadata),10,replace = F),]
    }
    
    # if(!all(grepl(".jpg", metadata$media_url, ignore.case=T))){
    #   cat("something's not an image", "\n")
    # }
    
      #--------------------------------------
      # Downloading images
      # Loop through the URLs and download the images
    # argument changes for download format
      for (i in seq_along(metadata$media_url)) {
        file_name <- paste0("virtual_herbarium/", 
                            paste0(gsub(" ","_",metadata$species[i]),"_",metadata$year[i],"_", metadata$key[i],".jpeg"))
        Sys.sleep(2)
        download.image(metadata$media_url[i], file_name)
        # resize?
        img <- image_read(file_name)
        image_write(img, file_name, quality = 10)
        cat("resized","\n")
      }
    }
}




  



# #----------------------------------
# {
#   # keep only specialist det
#   metadata <- subset(metadata, metadata$identifiedBy %in% specialists)
#   
#   # keep only those with year of collection in the metadata
#   metadata <- subset(metadata, !is.na(metadata$year))