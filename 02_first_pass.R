# Checking which species have at least 10 imaged specimens with valid coordinates on GBIF
# rm(list=ls())
# setwd("~/leaf_computer_vision")
library(mvh)

#-----------------------------
# Loading list of woody dicots from previous script.
woody_species <- read.csv("supporting_datasets/woody_species.csv")
woody_species <- woody_species$taxon_name

#-----------------------------
# Now searching to see which species have at least 10 specimens with valid coordinates on GBIF
list_metadata <- list()
for(i in 1:length(woody_species)) {
  list_metadata[[i]] <- search_specimen_metadata(taxon_name=woody_species[i], limit=10, hasCoordinate=T)
  cat(i, "\r")
}
save(list_metadata, file="supporting_datasets/list_metadata_complete.Rsave")
## NOTE: this file is too large to be committed and is on gitignore for now.

# load("supporting_datasets/list_metadata_complete.Rsave")
nrows <- c()
taxa_to_analyze <- c()
for(i in 1:length(list_metadata)) {
  nrows[i] <- nrow(list_metadata[[i]])
  cat(i, "\r")
  if(nrows[i]>=10) {
    taxa_to_analyze <- c(taxa_to_analyze, list_metadata[[i]]$scientificName[1])
  }
}

#-----------------------------
# This is the target list of taxa
write.csv(taxa_to_analyze, file="supporting_datasets/taxa_to_analyze.csv", row.names = F)

# filtered_list <- list_metadata
# filtered_list[nrows<10] <- NULL
# save(filtered_list, file="supporting_datasets/filtered_list.Rsave") # saving the original metadata too just in case

