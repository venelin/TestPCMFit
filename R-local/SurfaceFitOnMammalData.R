library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(TestPCMFit)
library(surface)


tree <- treeMammals
tree$edge.length <- tree$edge.length / max(PCMTreeNodeTimes(tree))
#tree$tip.label <- paste0("t", tree$tip.label)
#tree$node.label <- NULL
#tree <- nameNodes(tree)
values <- valuesMammals[, 1:PCMTreeNumTips(tree)]
values <- data.frame(t(values))
rownames(values) <- tree$tip.label

fitSurface <- runSurface(tree, values, plotaic = TRUE, verbose = TRUE)

#save(fitMappings, file = paste0("fitMappings_", id, "_t1_.RData"))


