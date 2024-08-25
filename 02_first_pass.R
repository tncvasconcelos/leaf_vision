setwd("/Users/tvasc/Desktop/leaf_computer_vision")
library(mvh)

woody_species <- read.csv("supporting_datasets/woody_species.csv")
woody_species <- woody_species$taxon_name

list_metadata <- list()
for(i in 1:length(woody_species)) {
  list_metadata[[i]] <- search_specimen_metadata(taxon_name=woody_species[i], limit=10, hasCoordinate=T)
  cat(i, "\r")
}
save(list_metadata, file="supporting_datasets/list_metadata_complete.Rsave")

#load("list_metadata_54915.Rsave")
nrows <- c()
taxa_to_analyze <- c()
for(i in 1:length(list_metadata)) {
  nrows[i] <- nrow(list_metadata[[i]])
  cat(i, "\r")
  if(nrows[i]>=10) {
    taxa_to_analyze <- c(taxa_to_analyze, list_metadata[[i]]$scientificName[1])
  }
}

write.csv(taxa_to_analyze, file="supporting_datasets/taxa_to_analyze.csv", row.names = F)

filtered_list <- list_metadata
filtered_list[nrows<10] <- NULL
save(filtered_list, file="supporting_datasets/filtered_list.Rsave")

