#LMA houwie
library(ape)
library(phylolm)
library(geiger)
library(phytools)
library(scales)
library(RColorBrewer)

setwd("~/leaf_vision/")

merged_dataset <- read.csv("data/merged_dataset.csv")
tre <- read.tree("trees/GBMB.tre")

lma_results <- aggregate(merged_dataset$LMA, list(merged_dataset$genus_species), 
  FUN = function(x) c(mean(log(x)), sd(log(x))/length(x)))
lma_results <- data.frame(sp = lma_results$Group.1,
  la = lma_results$x[,1],
  se =  lma_results$x[,2],
  row.names = lma_results$Group.1)


merged_dataset_2 <- merged_dataset[!duplicated(merged_dataset$filename),]
super_biomes <- data.frame(sp = merged_dataset_2$genus_species, 
  biome = as.factor(merged_dataset_2$super_biome))
super_biomes <- super_biomes[!is.na(super_biomes$biome),]
tmp <- aggregate(super_biomes$biome, by = list(super_biomes$sp), 
  FUN = function(x) as.numeric(table(factor(x, levels = levels(super_biomes$biome))) > 0))

# library(corHMM)
# 
# dat <- as.data.frame(do.call(cbind, tmp))
# colnames(dat) <- c("sp", levels(super_biomes$biome))
# missing_sp <- dat$sp[!dat$sp %in% tre$tip.label]
# dat <- dat[dat$sp %in% tre$tip.label, ]
# phy <- keep.tip(tre, dat$sp)
# phy <- ladderize(phy)
# phy$node.label <- NULL
# rownames(dat) <- dat$sp
# dat <- dat[phy$tip.label, ]
# lma_results <- lma_results[phy$tip.label, ]
# 
# lma_by_biome <- as.data.frame(aggregate(lma_results$la, 
#   by = list(dat$arid, dat$temperate, dat$tropical), mean))
# colnames(lma_by_biome) <- c(levels(super_biomes$biome), "lma")
# 
# dredge_fit <- corHMMDredge(phy, biome_df, 1, pen.type = "l1", lambda = 1)
# plotRECON(phy, test$states, show.tip.label = FALSE, type="fan")