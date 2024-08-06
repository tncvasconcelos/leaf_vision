# Downloading images and records from GBIF
# setwd("/Users/tvasc/Desktop/leaf_computer_vision")
source("00_functions.R")

# # Produce list of species to search:
# #-----------------------------
# # Load WCVP dataset
# dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
# 
# #-----------------------------
# # Merge them in one big table
# all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
# all_vars <- subset(all_vars, all_vars$taxon_rank=="Species")
# all_vars <- subset(all_vars, all_vars$taxon_status=="Accepted")
# unique_species <- unique(all_vars$taxon_name)
# write.csv(unique_species, file="species_to_search.csv", row.names = F)

# Let's look at the metadata first
unique_species <- read.csv("species_to_search.csv")
# ----------------------------------



all_names <- pbapply::pblapply(unique_species, gnr_resolve_x, cl=6)


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

#----------------------------------
full.herb.search_par <- function(names_to_search) {
  name_search_x <- function(x) {
    metadata <- NULL
    try(metadata <- full.herb.search(x), silent=T)
    if(is.null(metadata)) {
      metadata <- paste0(x,"_download_fail")
    }
    write.csv(metadata, file=paste0("metadata/", x, "_mvh_metadata.csv"), row.names=F)
    return(metadata)
  }
  all_names <- pbapply::pblapply(names_to_search, name_search_x, cl=6)
  return(as.list(all_names))
}

test<-full.herb.search_par(unique_species$x[1:1000])


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



