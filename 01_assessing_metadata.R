# Downloading images and records from GBIF
# rm(list=ls())
# setwd("/Users/tvasc/Desktop/leaf_computer_vision")
source("00_functions.R")

#-----------------------------
# Let's look at the metadata first
woody_species <- read.csv("supporting_datasets/woody_species.csv")
# ----------------------------------

sample_species <- sample(woody_species$taxon_name, 100)

name_search_x <- function(x) {
  metadata <- tryCatch({full.herb.search.metadata(x)}, error = function(e) {return(e$message)})
  write.csv(metadata, file=paste0("metadata/", x, "_mvh_metadata.csv"), row.names=F)
  return(metadata)
}

metadata <- pbapply::pblapply(sample_species, name_search_x)

all_media_links <- matrix(nrow=0, ncol=3)
for(i in 1:length(metadata)) {
  one_result <- metadata[[i]]
  if(class(one_result)=="data.frame"){
    all_media_links <- rbind(all_media_links, one_result[,c("species","key","media_url")])
  }
}

n_specimens <- 150
while(length(list.files("virtual_herbarium")) < n_specimens) {
  image_to_sample <- sample(1:nrow(all_media_links), 1)
  file_name <- paste0("virtual_herbarium/", 
                      paste0(gsub(" ","_",all_media_links$species[image_to_sample]),"_", all_media_links$key[image_to_sample],".jpeg"))
  Sys.sleep(2)
  try(download.herb.image(all_media_links$media_url[image_to_sample], file_name))
  try(try_img <- resize.image(file_name))
  if(exists("try_img")) {
    image_write(try_img, file_name)
    cat("resized","\n")
    remove("try_img")
  }
}


# Another way of resizing them:
# for(i in 1:length(list.files("virtual_herbarium"))) {
#   file_name <- list.files("virtual_herbarium", full.names = T)[i]
#   try(try_img <- image_read(file_name))
#   if(exists("try_img")) {
#     image_write(try_img, file_name, quality = 75)
#     cat("resized","\n")
#     remove("try_img")
#   }
# }




# For a full download:
# download_fail <- c()
# for(species_index in 1:length(unique_species)) {
#   metadata <- NULL
#   Sys.sleep(2)
#   try(metadata <- full.herb.search(unique_species[species_index]), silent=T)
#   if(!is.null(metadata)) {
#     #--------------------------------------
#     write.csv(metadata, file=paste0("metadata/",reference_table$gbif_name[species_index],".csv"), row.names=F)
#   } else {
#     download_fail <- c(download_fail, unique_species[species_index])
#     write.csv(download_fail, file="download_fail.csv", row.names=F)
#   }
# }





# de-duplicating a little with available information
# based on collector number
#subset_rec_number <- subset(metadata, !is.na(metadata$recordNumber))
#metadata <- subset(metadata, !metadata$key %in% subset_rec_number$key[duplicated(subset_rec_number$recordNumber[])])

# Filtering dataset a little before downloading
# remove things with no coordinates
metadata <- subset(metadata, !is.na(metadata$decimalLatitude))
metadata <- subset(metadata, !is.na(metadata$decimalLongitude))

metadata <- FilterWCVP_genus(points=metadata, all_vars, twgd_data, lon="decimalLongitude", lat="decimalLatitude") 


# # Load metadata
# all_metadata <- list.files("metadata", full.names = T)
# all_metadata <- lapply(all_metadata, read.csv)

#----------------------------------
# For just a look at the number of specimens available online
available_data <- matrix(nrow=nrow(reference_table), ncol=2)

for(species_index in 1:length(reference_table$wcvp_name)) { 
  metadata <- NULL
  try(metadata <- full.herb.search(reference_table$wcvp_name[species_index]), silent=T)
  if(!is.null(metadata)) {
    #--------------------------------------
    # Filtering dataset little before downloading
    # remove things with no coordinates
    metadata <- subset(metadata, !is.na(metadata$decimalLatitude))
    metadata <- subset(metadata, !is.na(metadata$decimalLongitude))
    
    # removing inaturalist observations
    metadata <- subset(metadata, metadata$basisOfRecord=="PRESERVED_SPECIMEN")
    metadata <- subset(metadata, !grepl("inaturalist",metadata$media_url))
    
    # removing records out of the native range
    #--------------------------------------
    #--------------------------------------
    available_data[species_index,1] <- reference_table$wcvp_name[species_index]
    available_data[species_index,2] <- nrow(metadata)
    write.csv(available_data, file="available_data.csv")
  } else {
    available_data[species_index,1] <- reference_table$wcvp_name[species_index]
    available_data[species_index,2] <- 0
  }
}
#----------------------------------



