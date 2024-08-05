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



