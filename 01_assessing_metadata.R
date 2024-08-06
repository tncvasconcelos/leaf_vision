# Downloading images and records from GBIF

setwd("/Users/tvasc/Desktop/leaf_computer_vision")
library(rgbif)
library(rvest)
library(httr)
library(dplyr)
library(tidyr)

source("00_functions.R")

trees <- load.trees("2_trees/")

# Pilot with Myrteae
tip_labels <- trees$Myrteae_NMWG_in_prep$tip.label # species sampled in tree
specialists <- c("Snow") # defining specialists?


#-----------------------------
# Load WCVP dataset
dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_vars <- subset(all_vars, all_vars$taxon_rank=="Species")
all_vars <- subset(all_vars, all_vars$taxon_status=="Accepted")
unique_species <- unique(all_vars$taxon_name)

reference_table <- list.files("taxized_reference_tables", full.names = T)
reference_table <- reference_table[grep("/reference_table",reference_table)]
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
reference_table <- subset(reference_table, reference_table$wcvp_name%in%unique_species)

# Let's look at the metadata first
#----------------------------------
# For a full download:
# for(species_index in 1:length(reference_table$gbif_name)) { 
#   metadata <- NULL
#   try(metadata <- full.herb.search(reference_table$gbif_name[species_index]), silent=T)
#   if(!is.null(metadata)) {
#     #--------------------------------------
#     # Filtering dataset little before downloading
#     # remove things with no coordinates
#     metadata <- subset(metadata, !is.na(metadata$decimalLatitude))
#     metadata <- subset(metadata, !is.na(metadata$decimalLongitude))
#     
#     # removing inaturalist observations
#     metadata <- subset(metadata, metadata$basisOfRecord=="PRESERVED_SPECIMEN")
#     metadata <- subset(metadata, !grepl("inaturalist",metadata$media_url))
#     
#     # de-duplicating a little with available information
#     # based on collector number
#     subset_rec_number <- subset(metadata, !is.na(metadata$recordNumber))
#     metadata <- subset(metadata, !metadata$key %in% subset_rec_number$key[duplicated(subset_rec_number$recordNumber[])])
#     
#     # reduce and standardize metadata
#     write.csv(metadata, file=paste0("metadata/",reference_table$gbif_name[species_index],".csv"), row.names=F)    
#   }
# }
#----------------------------------

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



