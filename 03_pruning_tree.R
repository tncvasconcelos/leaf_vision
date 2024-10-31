# Pruning Smith and Brown (2018) molecular tree to include only species that have at least 10 imaged specimens with valid coordinates
# rm(list=ls())
# setwd("/Users/tvasc/Desktop/leaf_computer_vision")

#-----------------------------
# Pruning tree to match list of species with 10 imaged specimens with valid coordinates
big_tree <- readRDS("supporting_datasets/taxized_GBMB.Rdata")
big_tree$tip.label <- unname(big_tree$tip.label)

# list of species with okay images and coordinates
species_to_sample <- read.csv("supporting_datasets/taxa_to_analyze.csv")
big_tree_pruned <- ape::keep.tip(big_tree, which(big_tree$tip.label %in% species_to_sample$x))
#ape::write.tree(big_tree_pruned, file="supporting_datasets/big_tree_pruned.tre")

# making a reference table because I messed up and wrote the taxized tree as newick
big_tree_pruned_newick <- phytools::read.newick("supporting_datasets/big_tree_pruned.tre")
big_tree_pruned_original <- big_tree_pruned

ref_table <- data.frame(taxized_names=big_tree_pruned_original$tip.label, newick_names=big_tree_pruned_newick$tip.label)
save(ref_table, file="supporting_datasets/ref_table.Rsave")

#-------
#big_tree$tip.label <- gsub(" .*","",big_tree$tip.label)
#genus_tree <- ape::drop.tip(big_tree, which(duplicated(big_tree$tip.label)))
