# Load the package and dependencies (ape, phytools, corpcor, subplex, spam)
library(PhylogeneticEM)
library(TestPCMFit)
library(PCMBase)

tree <- treeMammals
tree$edge.length <- tree$edge.length/max(PCMTreeNodeTimes(tree, tipsOnly = TRUE))
values <- valuesMammals


fitPhyloEM_scOU <- PhyloEM(tree, values, "scOU")

save(fitPhyloEM_scOU, file = "fitPhyloEM_scOU.RData")
