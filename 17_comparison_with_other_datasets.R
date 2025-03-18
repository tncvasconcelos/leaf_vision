##########################
# Revision
##########################
library(data.table)
library(raster)
library(BIEN)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

###################################
organize.bubble.plot <- function(trait_table, all_vars, twgd_data) {
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% trait_table[,1])
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  focal_areas <- unique(wcvp_subset$area_code_l3)
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

#########################
sum.twgd.trait <- function(one_dataset,twgd_data,all_vars) {
  tmp_rasters <- list()
  for(i in 1:nrow(one_dataset)) {
    wcvp_subset <- subset(all_vars, all_vars$taxon_name == one_dataset$species[i])
    occ_areas <- wcvp_subset$area_code_l3
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
    tmp0 <- raster(area_plus_buffer, res=0.1)
    tmp0[]  <- 1  
    tmp1 <- crop(tmp0, extent(area_plus_buffer))
    tmp2 <- mask(tmp1, area_plus_buffer)
    tmp_rasters[[i]] <- tmp2
    cat(paste0(i, " out of ", nrow(one_dataset)), "\r")
  }
  return(tmp_rasters)
}

#######
map.for.synthesis <- function(trait_table, all_vars, twgd_data) {
  trait_table=trait_table
  colnames(trait_table) <- c("species","gbif")
  organized_table_for_plot_total <- organize.bubble.plot(trait_table, all_vars, twgd_data)
  organized_table_for_plot_no_data <- organize.bubble.plot(subset(trait_table, trait_table$gbif=="no_data"), all_vars, twgd_data)
  
  nas <- rep(NA, nrow(organized_table_for_plot_no_data))
  proportion_table <- data.frame(sp_rich_prop=nas, sp_rich_total=nas, one_area=nas, lon=nas, lat=nas)
  for(i in 1:nrow(organized_table_for_plot_no_data)) {
    tmp_sp_rich <- organized_table_for_plot_no_data$sp_rich[i]
    one_area <- organized_table_for_plot_no_data$one_area[i]
    total_sp_rich <- organized_table_for_plot_total$sp_rich[organized_table_for_plot_total$one_area == one_area]
    one_proportion <- round(tmp_sp_rich / total_sp_rich, 3)
    proportion_table$sp_rich_prop[i] <- one_proportion
    proportion_table$sp_rich_total[i] <- total_sp_rich
    proportion_table$one_area[i] <- one_area
    proportion_table$lon[i] <- organized_table_for_plot_no_data$lon[i]
    proportion_table$lat[i] <- organized_table_for_plot_no_data$lat[i]
  }
  
  # hist(proportion_table$sp_rich_prop)
  # cutting the data
  
  return(proportion_table)
}

woody_species <- read.csv("supporting_datasets/woody_species_full_table.csv")

# making a "trait" dataset for species with no data in our dataset
merged_dataset <- read.csv("data/merged_dataset_final.csv")
lam_lm2_total <- gsub("_"," ",unique(merged_dataset$fullname))

lam_lm2_total <- setdiff(lam_lm2_total, lam_try_total)

lm2_gap_table <- data.frame(species=woody_species$taxon_name, lam_lm2=NA)
lm2_gap_table$lam_lm2[which(lm2_gap_table$species %in% lam_lm2_total)] <- "has_data"
lm2_gap_table$lam_lm2[which(!lm2_gap_table$species %in% lam_lm2_total)] <- "no_data"


organized_table_for_plot_total <- organize.bubble.plot(subset(lm2_gap_table, lm2_gap_table$lam_lm2=="has_data"), all_vars, twgd_data)
twgd_data_A <- merge(twgd_data01, organized_table_for_plot_total, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_A <- subset(twgd_data_A , twgd_data_A$LEVEL3_NAM!="Antarctica")
figure2_map_A <- ggplot(data = twgd_data_A) +
  geom_sf(aes(fill = sp_rich)) +
  scale_fill_viridis_c(option = "A") +
  theme_classic()
figure2_map_A



# Loading TRY data
try <- fread("R1/TRY/39869.txt")
#try1 <- head(try)
lam_try_total <- unique(try$SpeciesName)

try_gap_table <- data.frame(species=woody_species$taxon_name, lam_try=NA)
try_gap_table$lam_try[which(try_gap_table$species %in% lam_try_total)] <- "has_data"
try_gap_table$lam_try[which(!try_gap_table$species %in% lam_try_total)] <- "no_data"


###########
# Loading TDWG shape files (see README in data/wgsrpd-master/ for more information)
#-----------------------------
path="wcvp/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
twgd_data01 <- sf::st_as_sf(twgd_data)

# Figure 2A: molecular data for phylogenetic reconstruction
trait_table_A <- lm2_gap_table
proportion_table_A <- map.for.synthesis(trait_table=trait_table_A, all_vars=woody_species, twgd_data)
twgd_data_A <- merge(twgd_data01, proportion_table_A, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_A <- subset(twgd_data_A , twgd_data_A$LEVEL3_NAM!="Antarctica")
figure2_map_A <- ggplot(data = twgd_data_A) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic()

trait_table_B <- try_gap_table
proportion_table_B <- map.for.synthesis(trait_table=trait_table_B, all_vars=woody_species, twgd_data)
twgd_data_B <- merge(twgd_data01, proportion_table_B, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_B <- subset(twgd_data_B , twgd_data_B$LEVEL3_NAM!="Antarctica")
figure2_map_B <- ggplot(data = twgd_data_B) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic()



##########################
# Distribution of species in both datasets

# Distribution of species exclusive to our dataset



#-------------------
# MAPS WITH GLOBAL DISTRIBUTION
world <- ne_countries(scale = "medium", returnclass = "sf")
#------------------
# Comparing with BIEN
library(BIEN)

lam_bien <- BIEN_trait_trait(trait="leaf area per leaf dry mass")

lam_bien <- try

colnames(lam_bien)
#-------------------
lam_bien$grid_lat <- floor(lam_bien$latitude / 2.5) * 2.5
lam_bien$grid_lon <- floor(lam_bien$longitude / 2.5) * 2.5

grid_density <- lam_bien %>%
  group_by(grid_lat, grid_lon) %>%
  summarize(count = n(), .groups = "drop")

# Plot the heatmap
n_points_map <- ggplot(data = world) +
  geom_sf(fill = "lightgrey", color = "white") +  # World map base
  geom_tile(data = grid_density, 
            aes(x = grid_lon, y = grid_lat, fill = log(count)), 
            width = 2.5, height = 2.5, alpha=0.8) +  # 5x5 grid cells
  scale_fill_viridis_c(option = "plasma", name = "Point Count") +  # Density color scale
  theme_minimal() +
  labs(title = "Point Density Heatmap by 2.5x2.5 Degree Grid",
       x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "white"))


