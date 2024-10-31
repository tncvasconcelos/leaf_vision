# Creating phylogenetically balanced "random" sample of species
# rm(list=ls())
# setwd("~/leaf_computer_vision")

library(ape)

phy <- phytools::read.newick("supporting_datasets/big_tree_pruned.tre")

#pseudoinverse()
# C <- vcv(phy)
# C_inv<- solve(C)
# I <- diag(dim(C)[1])
# weights <- (t(I) %*% C_inv) %*% matrix(1, dim(C)[1], 1)
# weights_2 <- c(weights/sum(weights))

# Save the objects into a list
# data_list <- list(C_inv = C_inv, weights=weights)

# weights <- setNames(weights[,1], colnames(C_inv))
# Save the list as an RDS file
# saveRDS(data_list, file = "intermediate_datasets/phylo_data.rds")
# saveRDS(weights, file = "intermediate_datasets/weights.rds")
weights <- readRDS("supporting_datasets/weights.rds")

prelim_taxa_sample <- sample(names(weights), size = 2000, replace = FALSE, prob = weights)
write.csv(prelim_taxa_sample, "final_taxa_sample.csv")

par(mar=c(.1,.1,.1,.1))
plot(phy, show.tip.label = FALSE, type="fan")
tiplabels(pch = 16, col = "red", tip = match(prelim_taxa_sample, phy$tip.label), offset = 5)

