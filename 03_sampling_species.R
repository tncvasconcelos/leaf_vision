
setwd("/Users/tvasc/Desktop/leaf_computer_vision")

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

#-----------------------------
all_available_data <- read.csv("available_data.csv")
all_available_data <- subset(all_available_data, !is.na(all_available_data$V2))
all_available_data <- subset(all_available_data, all_available_data$V2 > 5)

#-----------------------------
# Finding out how much of each area has photos of specimens available
all_areas <- unique(all_vars$area_code_l3)
prop_sample <- as.data.frame(matrix(nrow=length(all_areas), ncol=3))
colnames(prop_sample) <- c("area_code","sp_rich","prop_available_sample")
for(area_index in 1:length(all_areas)) {
  prop_sample[area_index,1] <- all_areas[area_index]
  sp_rich0 <- all_vars$taxon_name[ all_vars$area_code_l3 == all_areas[area_index]]
  prop_sample[area_index,2] <- length(sp_rich0)
  sp_available <- sp_rich0[which(sp_rich0 %in% all_available_data$V1)]
  prop_sample[area_index,3] <- round((length(sp_available)) / (length(sp_rich0)),3)
}



