##########################
# Revision
# During revisions, reviewers asked to highlight how this approach can improve tropical sampling in comparison with what's already existing in datasets like TRY
##########################
library(data.table)
library(raster)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

###################################
organize.bubble.plot <- function(trait_table, all_vars, twgd_data) {
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% trait_table[,1])
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  twgd_data01 <- sf::st_as_sf(twgd_data)
  focal_areas <- as.character(unique(twgd_data01$LEVEL3_COD))
  results <- matrix(nrow=0, ncol=5)
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    one_subset <- subset(wcvp_subset, wcvp_subset$area_code_l3==one_area)
    sp_rich <- length(unique(one_subset$taxon_name))
    family_rich <- length(unique(one_subset$family))
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    if(nrow(area_plus_buffer)>0) {
      centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
      lon <- extent(centroids)[1]
      lat <- extent(centroids)[3]
      results <- rbind(results, cbind(sp_rich, family_rich, one_area, lon, lat))
    } 
    cat(i, "\r")
  }
  results <- as.data.frame(results)
  results$sp_rich <- as.numeric(results$sp_rich)
  results$family_rich <- as.numeric(results$family_rich)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}

# Load data about woody species from WCVP
woody_species <- read.csv("supporting_datasets/woody_species_full_table.csv")

###########
# Load TDWG shape files (see README in data/wgsrpd-master/ for more information)
#----------------------------- 
path="wcvp/wgsrpd-master/level3/level3.shp"
twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
twgd_data01 <- sf::st_as_sf(twgd_data)
#-----------------------------

# Loading TRY data
try <- fread("R1/TRY/39869.txt")
lam_try_total <- unique(try$SpeciesName)

# making a "trait" dataset for species with no data in our dataset
merged_dataset <- read.csv("data/merged_dataset_final.csv")
lam_lm2_total <- gsub("_"," ",unique(merged_dataset$fullname))

lam_lm2_total <- setdiff(lam_lm2_total, lam_try_total)

lm2_gap_table <- data.frame(species=woody_species$taxon_name, lam_lm2=NA)
lm2_gap_table$lam_lm2[which(lm2_gap_table$species %in% lam_lm2_total)] <- "has_data"
lm2_gap_table$lam_lm2[which(!lm2_gap_table$species %in% lam_lm2_total)] <- "no_data"

organized_table_for_plot_total <- organize.bubble.plot(trait_table=subset(lm2_gap_table, lm2_gap_table$lam_lm2=="has_data"), all_vars=woody_species, twgd_data)
twgd_data_A <- merge(twgd_data01, organized_table_for_plot_total, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_A <- subset(twgd_data_A , twgd_data_A$LEVEL3_NAM!="Antarctica")


library(dplyr)

# Replace 0 values with NA
twgd_data_A <- twgd_data_A %>%
  mutate(sp_rich = ifelse(sp_rich < 1, NA, sp_rich))

# Plot with white areas for 0 values
figure2_map_A <- ggplot(data = twgd_data_A) +
  geom_sf(aes(fill = sp_rich)) +
  scale_fill_viridis_c(
    option = "B", 
    na.value = "white", 
    direction = -1,               # Reverse the scale
    limits = c(min(twgd_data_A$sp_rich, na.rm = TRUE), 1.4 * max(twgd_data_A$sp_rich, na.rm = TRUE)),  # Exclude extremes
    #oob = scales::squish           # Squish out-of-bound values into the defined range    
    ) +
  theme_classic()

figure2_map_A



