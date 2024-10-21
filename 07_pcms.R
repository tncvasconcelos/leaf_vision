library(ape)
library(phylolm)
library(geiger)
library(phytools)
library(scales)
library(RColorBrewer)


setwd("~/leaf_vision/")

dat <- read.csv("data/lma_results.csv")
dat$specimen <- gsub("_\\d+$", "", dat$specimen)
dat$specimen <- gsub("_\\d+$", "", dat$specimen) # some had 2 numbers??
dat <- aggregate(dat$mean_LMA, by = list(dat$specimen), function(x) mean(log(x)))
plot(density(dat$x))
abline(v = mean(dat$x) - (3 * sd(dat$x)), col = "red", lwd = 2)
abline(v = mean(dat$x) + (3 * sd(dat$x)), col = "red", lwd = 2)

tre <- read.tree("trees/GBMB.tre")

missing_sp <- dat$Group.1[!dat$Group.1 %in% tre$tip.label]
dat <- dat[dat$Group.1 %in% tre$tip.label, ]

phy <- keep.tip(tre, dat$Group.1)
phy <- ladderize(phy)
rownames(dat) <- dat$Group.1
dat <- dat[phy$tip.label,]
H <- max(node.depth.edgelength(phy))

rescale_values <- function(x, a, b) {
  (a + (x - min(x)) / (max(x) - min(x)) * (b - a))
}
plot_tip_values <- setNames(rescale_values(dat$x, H*1.01, H*1.1), dat$Group.1)


lm_dat <- setNames(dat$x, dat$Group.1)
models <- c("BM", "lambda", "kappa", "EB")

all_fits <- lapply(models, function(x) phylolm(lm_dat~1, phy = phy, model = x, boot=100))
summary(all_fits[[2]])

myPalette <- brewer.pal("YlGnBu", n = 10)
asr <- fastAnc(phy, lm_dat)
node_vals <- c(setNames(dat$x, 1:Ntip(phy)), asr)
edge_vals <- c()
for(i in 1:nrow(phy$edge)){
  focal_edge <- phy$edge[i,]
  edge_vals <- c(edge_vals, mean(node_vals[focal_edge]))
}

plot.phylo(phy, show.tip.label = FALSE, no.margin = TRUE, direction = "upwards", y.lim=c(0, H*1.1))
for(i in seq_along(plot_tip_values)){
  segments(i, H*1.01, i, plot_tip_values[i])
}






