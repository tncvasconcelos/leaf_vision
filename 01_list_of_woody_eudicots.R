# Downloading images and records from GBIF
# setwd("/Users/tvasc/Desktop/leaf_computer_vision")
source("00_functions.R")

families_to_exclude <- read.csv("families_to_exclude.txt")
# Produce list of species to search:
#-----------------------------
# Load WCVP dataset
dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_vars <- subset(all_vars, all_vars$taxon_rank=="Species")
all_vars <- subset(all_vars, all_vars$taxon_status=="Accepted")
all_vars <- subset(all_vars, !all_vars$family%in%families_to_exclude$families_to_exclude)

#-----------------------------
# Filter to get just woody eudicots
life_form_scoring <- read.csv("lifeform_mapping.csv")
life_form_scoring <- rbind(life_form_scoring, c("",""))
life_forms <- subset(all_vars, !duplicated(all_vars$taxon_name))
life_forms$life_form <- NA 
for(i in 1:length(life_forms$lifeform_description)) {
  life_forms$life_form[i] <- life_form_scoring$humphreys_lifeform[which(life_form_scoring$lifeform_description==life_forms$lifeform_description[i])] 
  cat(i, "\r")
}
woody_species <- subset(life_forms, life_forms$life_form=="woody perennial")
unique_species <- unique(woody_species$taxon_name)
write.csv(unique_species, file="species_to_search.csv", row.names = F)
