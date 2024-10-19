library(ape)
library(phylolm)
library(geiger)
library(phytools)
library(scales)
library(RColorBrewer)

setwd("~/leaf_vision/")

dat <- read.csv("data/lma_results.csv")
tre <- read.tree("trees/GBMB.tre")

missing_sp <- dat$sp[!dat$sp %in% tre$tip.label]
dat <- dat[dat$sp %in% tre$tip.label, ]

phy <- keep.tip(tre, dat$sp)
phy <- ladderize(phy)
rownames(dat) <- dat$sp
dat <- dat[phy$tip.label,]
H <- max(node.depth.edgelength(phy))

rescale_values <- function(x, a, b) {
  (a + (x - min(x)) / (max(x) - min(x)) * (b - a))
}
plot_tip_values <- setNames(rescale_values(dat$lma, H*1.01, H*1.1), dat$sp)


lm_dat <- setNames(dat$lma, dat$sp)
models <- c("BM", "lambda", "kappa", "EB")

all_fits <- lapply(models, function(x) phylolm(lm_dat~1, phy = phy, measurement_error = dat$se, model = x, boot=1))
phylolm(lm_dat~1, phy  = phy, measurement_error = TRUE, model = "BM")
summary(all_fits[[2]])

myPalette <- brewer.pal("YlGnBu", n = 9)
asr <- fastAnc(phy, lm_dat)
all_nodes <- c(setNames(dat$lma, 1:Ntip(phy)), asr)
col_func <- colorRampPalette(c("darkblue", "darkred"))
cols <- col_func(100)
bins <- cut(all_nodes, breaks = 100)
all_cols <- cols[match(bins, levels(bins))]

plot.phylo(phy, show.tip.label = FALSE, no.margin = TRUE, direction = "upwards", 
  y.lim=c(0, H*1.1), edge.width = 0.5, edge.color = all_cols[phy$edge[,1]])
for(i in seq_along(plot_tip_values)){
  segments(i, H*1.01, i, plot_tip_values[i], col = all_cols[i], lwd = 0.5)
}

