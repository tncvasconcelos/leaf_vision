library(taxize)
library(ape)

resolve.names <- function(names_to_solve) {
  gnr_resolve_x <- function(x) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=TRUE)$matched_name)
    if(is.null(tmp.name)) {
      tmp.name <- paste0(x,"_UNMATCHED")
    }
    return(tmp.name)
  }
  all_names <- pbapply::pblapply(names_to_solve, gnr_resolve_x, cl=6)
  return(as.character(all_names))
}

load.trees <- function(tree.dir) {
  tree_files <- list.files(tree.dir, full.names = T)
  all_trees <- list()
  for(i in 1:length(tree_files)) {
    load(tree_files[i])
    if(exists("one_tree")) {
      all_trees[[i]] <- one_tree
      names(all_trees)[i] <- gsub(paste0(c(paste0(tree.dir,"/"), ".Rsave"), collapse="|"),"", tree_files[i])
      rm("one_tree")
    }
  }
  return(all_trees)
}

#
full.herb.search <- function(species_name) {
  #--------------------------------------
  # Search GBIF for records with images
  all_gbif_data <- occ_search(scientificName = species_name , mediaType = "StillImage")
  
  #--------------------------------------
  # Extract URL and licence types
  metadata <- as.data.frame(all_gbif_data$data)
  metadata$media_url <- NA
  metadata$license <- NA
  for(obs_index in 1:nrow(metadata)) {
    media_info <- all_gbif_data$media[[obs_index]][[1]]
    metadata$media_url[obs_index] <- media_info[[1]]$identifier
    metadata$license <- media_info[[1]]$license
  }
  cat("Search for", species_name, "done!", "\n")
  return(metadata)
}

# Function to download images
download.image <- function(url, destfile) {
  tryCatch({
    download.file(url, destfile, mode = "wb")
  }, error = function(e) {
    message("Error downloading ", url)
  })
}


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

