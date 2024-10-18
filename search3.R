#---------------
# Downloading herbarium specimens
setwd("/Users/tvasc/Desktop/leaf_computer_vision")
# rm(list=ls())
library(mvh)

species_to_sample <- read.csv("final_taxa_sample.csv")[,2]
load("supporting_datasets/ref_table.Rsave")
taxized_names <- unique(ref_table[,1][ref_table$newick_names %in% species_to_sample])


#first_batch <- taxized_names[1:100]
for(i in 1:length(taxized_names)) {
  metadata <- search_specimen_metadata(taxon_name=taxized_names[i], limit=100, hasCoordinate=T)
  # Removing duplicates based on coordinates
  metadata <- subset(metadata, !duplicated(paste0(metadata$decimalLatitude,"_",metadata$decimalLongitude)))
  if(nrow(metadata)>5) {
    # Adding preference scores
    metadata$preference_score <- 0
    metadata$preference_score <-  metadata$preference_score + metadata$year
    # Favor eventDate not between September and May for decimalLatitude > 30
    metadata$preference_score <- metadata$preference_score + ifelse(metadata$decimalLatitude > 30 & (metadata$month >= 9 | metadata$month <= 5), 0, 1000)
    # Favor eventDate between September and May for decimalLatitude < 30
    metadata$preference_score <- metadata$preference_score + ifelse(metadata$decimalLatitude < 30 & (metadata$month >= 9 & metadata$month <= 5), 1000, 0)
    # A lot of these are field photos, so let's reduce their score
    metadata$preference_score <- metadata$preference_score - ifelse(grepl("urn:catalog:MO:Tropicos", metadata$identifier), 1000, 0)
    if(!is.null(metadata$institutionCode)) {
      # NYBG scores higher
      metadata$preference_score <- metadata$preference_score + ifelse(metadata$institutionCode%in%"NY", 1000, 0)
      # BRI,NSW and RBGE score lower (low res)
      metadata$preference_score <- metadata$preference_score - ifelse(metadata$institutionCode%in%"NSW", 1000, 0)
      metadata$preference_score <- metadata$preference_score - ifelse(metadata$institutionCode%in%"BRI", 1000, 0)
      metadata$preference_score <- metadata$preference_score - ifelse(metadata$institutionCode%in%"RBGE", 1000, 0)
      metadata$preference_score <- metadata$preference_score - ifelse(metadata$institutionCode%in%"P", 1000, 0)
    }
    # Sort by preference score (higher is better)
    metadata <- metadata[order(metadata$preference_score,decreasing=T), ]
    # Sample up to 5 rows from the sorted data.frame
    metadata <- head(metadata, 20)
  }
  write.csv(metadata, file=paste0("virtual_herbarium_NPleafPaper/metadata/",metadata$scientificName[1]," metadata.csv"))
  
  n_images_in_folder <- 0
    for(j in 1:nrow(metadata)) {  
      try(download_specimen_images(metadata[j,],
                                   dir_name="virtual_herbarium_NPleafPaper/final_dataset/",
                                   result_file_name=paste0("virtual_herbarium_NPleafPaper/result_files/",metadata$scientificName[j],metadata$gbifID[j])))
      
      # check quality
      file_name <- paste0("virtual_herbarium_NPleafPaper/final_dataset/",gsub(" ","_",paste(metadata$species[j],metadata$gbifID[j])),".jpeg")
      try_img <- try(magick::image_read(file_name), silent = TRUE)
      if(!inherits(try_img, "try-error")) {  # Check if resizing succeeded
        current_width <- magick::image_info(try_img)$width
        current_height <- magick::image_info(try_img)$height
        # Calculate the current megapixels
        current_megapixels <- (current_width * current_height) / 1e6
        if(length(current_megapixels) == 0) {
          file.remove(file_name)
          next
        }
        if(current_megapixels < 3) {
          file.remove(file_name)
        }
      } else {
        file.remove(file_name)
      }
      n_images_in_folder <- length(list.files("virtual_herbarium_NPleafPaper/final_dataset/",gsub(" ","_",paste(metadata$species[1]))))
      if(n_images_in_folder>=5){
        break
      }
    }   
  cat(i, "\r")
}

