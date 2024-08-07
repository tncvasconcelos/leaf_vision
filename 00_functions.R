library(taxize)
library(ape)
library(rgbif)
library(rvest)
library(httr)
library(dplyr)
library(tidyr)
library(magick)

########################################
# mvh functions
########################################

full.herb.search.metadata <- function(species_name) {
  #--------------------------------------
  # Search GBIF for records with images
  Sys.sleep(2)
  all_gbif_data <- occ_search(scientificName = species_name , mediaType = "StillImage")
  #--------------------------------------
  # Extract URL and licence types
  metadata <- as.data.frame(all_gbif_data$data)
  metadata$media_url <- NA
  metadata$license <- NA
  for(obs_index in 1:nrow(metadata)) {
    media_info <- all_gbif_data$media[[obs_index]][[1]]
    if("identifier" %in% names(media_info[[1]])) {
      metadata$media_url[obs_index] <- media_info[[1]]$identifier
      metadata$license <- media_info[[1]]$license
    }
  }
  #cat("Search for", species_name, "done!", "\n")
  metadata <- subset(metadata, metadata$basisOfRecord=="PRESERVED_SPECIMEN")
  metadata <- subset(metadata, !grepl("inaturalist",metadata$media_url))
  return(metadata)
}

# Function to download images
download.herb.image <- function(url, destfile) {
  tryCatch({
    download.file(url, destfile, mode = "wb")
  }, error = function(e) {
    message("Error downloading ", url)
  })
}

#----------------------------------
full.herb.search <- function(names_to_search, download_metadata=T, download_specimens=T, resize=T, n_cores=6) {
    metadata <- tryCatch({full.herb.search.metadata(names_to_search)}, error = function(e) {return(e$message)})
    if(download_metadata) {
      # cdd code to create a medatata folder first
      write.csv(metadata, file=paste0("metadata/", names_to_search, "_mvh_metadata.csv"), row.names=F)
    }
    if(download_specimens) {
      for (i in seq_along(metadata$media_url)) {
        file_name <- paste0("virtual_herbarium/", 
                            paste0(gsub(" ","_",metadata$species[i]),"_", metadata$key[i],".jpeg"))
        Sys.sleep(2)
        download.herb.image(metadata$media_url[i], file_name)
        # resize?
        if(resize) {
          try(try_img <- image_read(file_name))
          if(exists("try_img")) {
            image_write(try_img, file_name, quality = 20)
            cat("resized","\n")
            remove("try_img")
          }
        } 
      }
    }
  return(metadata)
}


########################################
# Functions to filter data
########################################
FilterWCVP_genus <- function(points, all_vars, twgd_data, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  #twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  
  wcvp_subset <- subset(all_vars, all_vars$genus %in% tmp_points$genus)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  occ_areas <- wcvp_subset$area_code_l3
  area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
  if(nrow(area_plus_buffer)>0) {
    coords <- tmp_points[,c("x","y")]
    sp::coordinates(coords) <- ~ x + y
    answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
    if(length(answer) != 0) {
      dubiousGBIF_ids <- as.character(tmp_points$gbifID[answer])
    }
  }  
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}


FilterWCVP <- function(points, all_vars, reference_table, twgd_data, species= "scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  #twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  for(species_index in 1:nrow(reference_table)) {
    gbif_subset <- subset(tmp_points, tmp_points$scientificName == reference_table$gbif_name[species_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
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
          dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
        }
      }
    }
    cat(species_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}

resize.image <- function(image_path, min_megapixels = 20, max_megapixels = 25) {
  # Load the image
  img <- image_read(image_path)
  
  # Get the current dimensions of the image
  current_width <- image_info(img)$width
  current_height <- image_info(img)$height
  
  # Calculate the current megapixels
  current_megapixels <- (current_width * current_height) / 1e6
  
  # Check if the current megapixels are already within the desired range
  if (current_megapixels >= min_megapixels && current_megapixels <= max_megapixels) {
    return(img)
  }
  
  # Calculate the scaling factor needed to get within the desired megapixel range
  scaling_factor <- sqrt(min_megapixels / current_megapixels)
  scaled_width <- round(current_width * scaling_factor)
  scaled_height <- round(current_height * scaling_factor)
  
  # Ensure the scaled image is not too large
  if (scaled_width * scaled_height / 1e6 > max_megapixels) {
    scaling_factor <- sqrt(max_megapixels / current_megapixels)
    scaled_width <- round(current_width * scaling_factor)
    scaled_height <- round(current_height * scaling_factor)
  }
  
  # Resize the image
  resized_img <- image_resize(img, geometry_size_pixels(scaled_width, scaled_height, preserve_aspect = TRUE))
  
  return(resized_img)
}
