#---------------
# Downloading herbarium specimens
setwd("/Users/tvasc/Desktop/leaf_computer_vision")
# rm(list=ls())
library(mvh)

species_to_sample <- read.csv("supporting_datasets/prelim_taxa_sample.csv")[,2]
load("supporting_datasets/ref_table.Rsave")
taxized_names <- unique(ref_table[,1][ref_table$newick_names %in% species_to_sample])


for(i in 1:length(taxized_names[1:100])) {
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
    # NYBG scores higher
    metadata$preference_score <- metadata$preference_score + ifelse(metadata$institutionCode%in%"NY", 1000, 0)
    # BRI and NSW score lower (low res)
    metadata$preference_score <- metadata$preference_score - ifelse(metadata$institutionCode%in%"NSW", 1000, 0)
    metadata$preference_score <- metadata$preference_score - ifelse(metadata$institutionCode%in%"BRI", 1000, 0)
    # Sort by preference score (higher is better)
    metadata <- metadata[order(metadata$preference_score,decreasing=T), ]
    # Sample up to 5 rows from the sorted data.frame
    metadata <- head(metadata, 5)
  }
  
  write.csv(metadata, file=paste0("virtual_herbarium_NPleafPaper/metadata/",metadata$scientificName[1]," metadata.csv"))
  try(download_specimen_images(metadata,
                           dir_name="virtual_herbarium_NPleafPaper/search1_100",
                           result_file_name=paste0("virtual_herbarium_NPleafPaper/search1_100/",metadata$scientificName[1])))
  
  cat(i, "\r")
}




